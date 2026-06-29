library(tidyverse)
library(randomForest)
library(quantreg)        # rq
library(grf)             # quantile_forest
library(caret)           # CV folds
library(SuperLearner)    # SuperLearner
library(tictoc)
library(glue)
library(future)
library(furrr)
library(ranger)
library(BART)
library(BAS)
RNGkind("L'Ecuyer-CMRG")

## ---- set working directory based on this script path ----
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args0[grep("^--file=", args0)])
if (length(script_path)) {
  setwd(dirname(normalizePath(script_path)))
}
# functions --------------------------------------------------------------------
source("0_functions_sim.R")

# global parameters
alpha_overall <- 0.15 # for once coverage 80% and for making decision on conformal p-values (same as gamma) 
alpha_FDR <- 0.15
SL_LIB_AR     <- c("SL.glm", "SL.ranger")
SL_LIB_NUIS   <- c("SL.glm", "SL.ranger")
seed <- 2025
n_rep <- 1000

w0 <- 0.14
b0 <- 0.01
n_max <- 100 # number of testing subjects 
decay <- 0.3
xi_long_mod <- generate_xi_mod(n_max, alpha = alpha_FDR, b0 = b0, decay=decay)

# for running each grid of scenario
run_grid <- function(
    cases        = c(1, 2),
    caps         = c(FALSE, TRUE),
    cs_types     = c("res", "CQR","clip"),
    # cvs          = c(FALSE, TRUE),
    cvs_nuisance = c(FALSE, TRUE),
    cvs_model    = c(FALSE, TRUE),
    # methods      = c("ID_reg","ID_rf","ID_sl","ED_reg","ED_rf","ED_sl"),
    methods = c("ID_reg","ID_rf","ID_sl","ID_bart","ID_bma","ED_reg","ED_rf","ED_sl","ED_bart","ED_bma"),
    gamma_thresh = alpha_overall,
    parallel     = FALSE,
    workers      = NULL,
    n_rep        = NULL,
    seed_rep     = 2025,
    strategy     = c("multisession","multicore")
){
  strategy <- match.arg(strategy)
  all_results <- list()
  summary_rows <- list()
  
  for (case in cases){
    # prepare data once per case
    prepared <- prepare_splits_for_case(case)
    split_ID <- prepared$split_ID
    split_ED <- prepared$split_ED
    sim_test <- prepared$sim_test
    
    if (!is.null(n_rep)) {
      total_rep <- length(split_ID)
      if (n_rep > total_rep) {
        stop("n_rep cannot be larger than total number of reps: ", total_rep)
      }
      set.seed(seed_rep)
      idx <- sort(sample.int(total_rep, n_rep))   # random but reproducible
      
      split_ID <- split_ID[idx]
      split_ED <- split_ED[idx]
      sim_test <- sim_test[idx]
    }
    
    for (cap in caps){
      for (cs_type in cs_types){
        for (cv_model in cvs_model) {
          for (cv_nuisance in cvs_nuisance){
            res <- run_methods_for_scenario(
              methods = methods,
              split_ID = split_ID, split_ED = split_ED, sim_test = sim_test,
              cs_type = cs_type, cv_nuisance = cv_nuisance, cv_model=cv_model, cap = cap,
              parallel = parallel, workers = workers, strategy = strategy
            )
            # summ <- summarize_methods(res, sim_test, gamma = gamma_thresh,
            #                           alpha_FDR = alpha_FDR, w0 = w0, b0 = b0, xi_mod = xi_long_mod)
            summ <- summarize_methods(
              results_list = res,
              sim_test     = sim_test,
              gamma        = gamma_thresh,
              alpha_FDR    = alpha_FDR,
              FDR_methods  = c("None","BH","Bonferroni","mLORD"),
              w0           = w0,
              b0           = b0,
              decay        = decay,
              xi_mod       = xi_long_mod
            )
            
            key <- glue("case{case}__cap{cap}__cs{cs_type}__cv_nuisance{cv_nuisance}__cv_model{cv_model}")
            all_results[[key]] <- res
            
            summary_rows[[key]] <- summ$wide %>%
              mutate(case = case, cap = cap, cs_type = cs_type, cv_nuisance = cv_nuisance, cv_model = cv_model, .before = 1)
          }
        }
      }
    }
  }
  
  summary_tidy <- bind_rows(summary_rows)
  # optional: also make a nice printed table per scenario
  list(
    raw = all_results,
    summary = summary_tidy %>% relocate(method, .after = last_col()),
    summary_pretty = summary_tidy %>%
      mutate(across(where(is.numeric) & !any_of(c("case","cap","cv_nuisance","cv_model")), ~round(.x, 4)))
  )
}
## ---- make results + logs directories ----
if (!dir.exists("results_cv_online")) dir.create("results_cv_online", showWarnings = FALSE)
if (!dir.exists("logs_cv_online"))    dir.create("logs_cv_online",    showWarnings = FALSE)

## ---- read Slurm array ID to choose scenario ----
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
cat("SLURM_ARRAY_TASK_ID =", task_id, "\n")

config <- expand.grid(
  case   = 1:36, # 1:6 (rct=50); 7:12 (rct=200)
  ctype  = "clip", # c("res", "clip")
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  arrange(case, ctype)
if (task_id < 1 || task_id > nrow(config)) {
  stop("task_id out of range. Must be between 1 and ", nrow(config), ".")
}
this_case  <- config$case[task_id]
this_ctype <- config$ctype[task_id]
scenario <- glue("case{this_case}_{this_ctype}")
cat("Running scenario:", scenario, "\n")

n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
# n_workers <- max(1, parallel::detectCores() - 2)
n_workers <- max(1, n_workers - 1)
n_workers <- min(n_workers, 6)
# plan(multicore, workers = n_workers)
future::plan(future::multisession, workers = n_workers)

grid_out <- run_grid(
  cases    = this_case,
  caps     = FALSE,         # keep whatever you used before
  cs_types = this_ctype,    # "clip"
  # cvs      = TRUE,
  cvs_model = TRUE, # previous FALSE
  cvs_nuisance = TRUE,
  parallel = TRUE,
  n_rep    = n_rep,
  seed_rep = seed,
  workers  = n_workers,
  strategy = "multisession"#"multicore"
)

rdata_file <- file.path("results_cv_online", glue("grid_out_{scenario}.RData"))
save(grid_out, file = rdata_file)

rdata_file_summary  <- file.path("results_cv_online",  glue("grid_out_{scenario}_summary.RData"))
summary_pretty <- grid_out$summary_pretty
save(summary_pretty, file = rdata_file_summary)
