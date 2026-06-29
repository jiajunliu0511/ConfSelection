library(tidyverse)
library(randomForest)
library(quantreg)        
library(grf)             
library(caret)           
library(SuperLearner)    
library(tictoc)
library(glue)
library(future)
library(furrr)
library(ranger)
library(BART)
library(BAS)
RNGkind("L'Ecuyer-CMRG")

args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args0[grep("^--file=", args0)])
if (length(script_path)) {
  setwd(dirname(normalizePath(script_path)))
}

source("0_functions_sim.R")

# global parameters
alpha_overall <- 0.15 
alpha_FDR <- 0.15
SL_LIB_AR <- c("SL.glm", "SL.ranger")
SL_LIB_NUIS <- c("SL.glm", "SL.ranger")
seed <- 2025

n_rep <- 200 # for CDF curves

# prepare data* retain Y1 and Y0 in calib
split_data <- function(data, a, ed, use_ed, borrow_trt, split_train = 0.5, split_nuisance = 0.2, method = "DR"){
  dat <- list(X=data$X[data$A==a,], Y=data$Y[data$A==a], Y1=data$Y1[data$A==a], Y0=data$Y0[data$A==a], S=data$S[data$A==a])
  id_rct <- seq_along(dat$Y)
  
  if (method == "DR"){
    n_nuisance <- ceiling(split_nuisance * length(dat$Y))
    n_train <- ceiling(split_train * (length(dat$Y) - n_nuisance))
    n_cal <- length(dat$Y) - n_train - n_nuisance
    split_id <- split(id_rct, sample(c(rep(1, n_train), rep(2, n_cal), rep(3, n_nuisance))))
    id_train <- split_id[[1]]; id_cal <- split_id[[2]]; id_nuisance <- split_id[[3]]
    
    dat_calib <- tibble(Y=dat$Y[id_cal], Y1=dat$Y1[id_cal], Y0=dat$Y0[id_cal], X=dat$X[id_cal,,drop=FALSE], S=dat$S[id_cal])
    dat_nuisance <- tibble(Y=dat$Y[id_nuisance], X=dat$X[id_nuisance,,drop=FALSE], S=dat$S[id_nuisance])
  } else {
    n_train <- ceiling(split_train * length(dat$Y))
    n_cal <- length(dat$Y) - n_train
    split_id <- split(id_rct, sample(c(rep(1, n_train), rep(2, n_cal))))
    id_train <- split_id[[1]]; id_cal <- split_id[[2]]

    dat_calib <- tibble(Y=dat$Y[id_cal], Y1=dat$Y1[id_cal], Y0=dat$Y0[id_cal], X=dat$X[id_cal,,drop=FALSE], S=dat$S[id_cal])
    dat_nuisance <- NULL
  }
  if (!use_ed){
    dat_train <- tibble(Y = dat$Y[id_train], X = dat$X[id_train, , drop = FALSE], S = dat$S[id_train])
  } else {
    if (a == 0){
      ed_train <- tibble(Y = ed$Y[ed$A == a], X = ed$X[ed$A == a, , drop = FALSE], S = ed$S[ed$A == a])
      rct_train <- tibble(Y = dat$Y[id_train], X = dat$X[id_train, , drop = FALSE], S = dat$S[id_train])
      dat_train <- bind_rows(rct_train, ed_train)
    } else {
      if (isTRUE(borrow_trt)){
        ed_train <- tibble(Y = ed$Y[ed$A == a], X = ed$X[ed$A == a, , drop = FALSE], S = ed$S[ed$A == a])
        rct_train <- tibble(Y = dat$Y[id_train], X = dat$X[id_train, , drop = FALSE], S = dat$S[id_train])
        dat_train <- bind_rows(rct_train, ed_train)
      } else {
        dat_train <- tibble(Y = dat$Y[id_train], X = dat$X[id_train, , drop = FALSE], S = dat$S[id_train])
      }
    }
  }
  lst(dat_train, dat_calib, dat_nuisance)
}


Extract_StochScores_DR <- function(split_data, ps = 0.5,
                                   fit_nuisance = "rf", fit_train = "rf", 
                                   CV_nusiance = TRUE, CV_model = TRUE, CV_fold = 5,
                                   conformal_type = "clip", cap = FALSE,
                                   SL_lib_nuisance = SL_LIB_NUIS, SL_lib_ar = SL_LIB_AR) {
  
  # 1. Bind Data (Strictly RCT data, A=1 and A=0)
  dat_nuisance <- bind_rows(split_data$split_dat1$dat_nuisance %>% mutate(A = 1),
                            split_data$split_dat0$dat_nuisance %>% mutate(A = 0)) %>% mutate(type = "nuisance")
  dat_calib <- bind_rows(split_data$split_dat1$dat_calib %>% mutate(A = 1),
                         split_data$split_dat0$dat_calib %>% mutate(A = 0)) %>% mutate(type = "calibration")
  dat_train <- bind_rows(split_data$split_dat1$dat_train %>% mutate(A = 1),
                         split_data$split_dat0$dat_train %>% mutate(A = 0)) %>% mutate(type = "train")
  
  dat_all <- bind_rows(dat_nuisance, dat_train, dat_calib)
  df_all <- df_X(dat_all)
  X_cols <- df_all %>% select(starts_with("X")) %>% names()
  
  ps_i <- rep(ps, nrow(df_all)) 
  
  # 2. Nuisance Modeling & Pseudo-Outcomes
  if (!CV_nusiance) {
    df_nuisance <- df_all %>% filter(type == "nuisance")
    train1 <- df_nuisance %>% filter(A == 1) %>% select(Y, all_of(X_cols))
    train0 <- df_nuisance %>% filter(A == 0) %>% select(Y, all_of(X_cols))
    
    hatmu1 <- fit_predict_gaussian(train1, df_all, "Y", X_cols, fit_nuisance, SL_lib_nuisance)
    hatmu0 <- fit_predict_gaussian(train0, df_all, "Y", X_cols, fit_nuisance, SL_lib_nuisance)
    
    A <- df_all$A
    Y <- df_all$Y
    df_all$Y_pseudo <- (A - ps_i) / (ps_i * (1 - ps_i)) * (Y - (hatmu0 * (1 - A) + hatmu1 * A)) + hatmu1 - hatmu0
    
    df_train <- df_all %>% filter(type == "train")
    df_calib_subset <- df_all %>% filter(type == "calibration")
  } else {
    df_model <- df_all %>% dplyr::filter(type %in% c("train","nuisance"))
    df_calib_subset <- df_all %>% dplyr::filter(type == "calibration")
    
    K <- CV_fold 
    folds <- caret::createFolds(df_model$Y, k = K)
    Y_pseudo_model <- rep(NA_real_, nrow(df_model))
    Y_pseudo_calib_mat <- matrix(NA_real_, nrow = nrow(df_calib_subset), ncol = K)
    
    for (k in seq_len(K)) {
      fold_idx  <- folds[[k]]
      train_idx <- setdiff(seq_len(nrow(df_model)), fold_idx)
      df_train_k <- df_model[train_idx, , drop = FALSE]
      
      train1_k <- df_train_k %>% dplyr::filter(A == 1) %>% dplyr::select(Y, dplyr::all_of(X_cols))
      train0_k <- df_train_k %>% dplyr::filter(A == 0) %>% dplyr::select(Y, dplyr::all_of(X_cols))
      
      df_fold_pred <- df_model[fold_idx, , drop = FALSE]
      hatmu1_k <- fit_predict_gaussian(train1_k, df_fold_pred, "Y", X_cols, fit_nuisance, SL_lib_nuisance)
      hatmu0_k <- fit_predict_gaussian(train0_k, df_fold_pred, "Y", X_cols, fit_nuisance, SL_lib_nuisance)
      
      hatmu1_calib <- fit_predict_gaussian(train1_k, df_calib_subset, "Y", X_cols, fit_nuisance, SL_lib_nuisance)
      hatmu0_calib <- fit_predict_gaussian(train0_k, df_calib_subset, "Y", X_cols, fit_nuisance, SL_lib_nuisance)
      
      A_k <- df_model$A[fold_idx]
      Y_k <- df_model$Y[fold_idx]
      Y_pseudo_model[fold_idx] <- (A_k - ps) / (ps * (1 - ps)) * (Y_k - (hatmu0_k * (1 - A_k) + hatmu1_k * A_k)) + (hatmu1_k - hatmu0_k)
      
      A_c <- df_calib_subset$A
      Y_c <- df_calib_subset$Y
      Y_pseudo_calib_mat[, k] <- (A_c - ps) / (ps * (1 - ps)) * (Y_c - (hatmu0_calib * (1 - A_c) + hatmu1_calib * A_c)) + (hatmu1_calib - hatmu0_calib)
    }
    df_train <- df_model
    df_train$Y_pseudo <- Y_pseudo_model
    df_calib_subset$Y_pseudo <- rowMeans(Y_pseudo_calib_mat)
  }
  
  # 3. Fit CATE Model
  if(!CV_model){
    mu_cal <- fit_predict_gaussian(df_train, df_calib_subset, "Y_pseudo", X_cols, fit_train, SL_lib_ar)
  } else { 
    K_pred <- CV_fold
    folds_model <- caret::createFolds(df_train$Y_pseudo, k = K_pred)
    mu_cal_mat <- matrix(NA_real_, nrow = nrow(df_calib_subset), ncol = K_pred)
    
    for (k in seq_len(K_pred)) {
      fold_idx_pred  <- folds_model[[k]]
      train_idx_pred <- setdiff(seq_len(nrow(df_train)), fold_idx_pred)
      df_train_k_pred <- df_train[train_idx_pred, , drop = FALSE]
      
      mu_cal_mat[, k] <- fit_predict_gaussian(df_train_k_pred, df_calib_subset, "Y_pseudo", X_cols, fit_train, SL_lib_ar)
    }
    mu_cal <- rowMeans(mu_cal_mat)
  }
  
  # ITEs
  True_ITE   <- df_calib_subset$Y1 - df_calib_subset$Y0
  Pseudo_ITE <- df_calib_subset$Y_pseudo
  Fitted_ITE <- mu_cal
  
  if (conformal_type == "res"){
    if(cap){
      Score_Oracle <- pmax(0, True_ITE - Fitted_ITE)
      Score_Pseudo <- pmax(0, Pseudo_ITE - Fitted_ITE)
    }else{
      Score_Oracle <- (True_ITE - Fitted_ITE)
      Score_Pseudo <- (Pseudo_ITE - Fitted_ITE)
    } 
  } else if (conformal_type == "clip"){
    M <- 100
    if(cap){
      Score_Oracle <- pmax(0, M * (True_ITE > 0) - Fitted_ITE)
      Score_Pseudo <- pmax(0, M * (Pseudo_ITE > 0) - Fitted_ITE)
    }else{
      Score_Oracle <- M * (True_ITE > 0) - Fitted_ITE
      Score_Pseudo <- M * (Pseudo_ITE > 0) - Fitted_ITE
    }   
  }
  
  tibble::tibble(
    Subject_ID   = 1:nrow(df_calib_subset),
    True_ITE     = True_ITE,
    Fitted_ITE   = Fitted_ITE,
    Pseudo_ITE   = Pseudo_ITE,
    Score_Oracle = Score_Oracle,
    Score_Pseudo = Score_Pseudo
  )
}

run_grid_stoch <- function(cases, ctype, cvs_nuisance=TRUE, cvs_model=TRUE, 
                           methods = c("ID_reg","ID_rf","ID_sl","ID_bart","ID_bma",
                                       "ED_reg","ED_rf","ED_sl","ED_bart","ED_bma"), workers) {
  
  all_scores <- list()
  
  for (case in cases) {
    prepared <- prepare_splits_for_case(case)
    split_ID <- prepared$split_ID[1:n_rep]
    split_ED <- prepared$split_ED[1:n_rep]
    sim_test <- prepared$sim_test[1:n_rep]
    
    for (m in methods) {
      if (grepl("^ED_", m)) {
        cfg <- method_mapper(m, split_ED, split_ID) # Train on ED, Calibrate on ID
      } else {
        cfg <- method_mapper(m, split_ID, split_ID) # Train on ID, Calibrate on ID
      }
      
      res_list <- future_map(cfg$split, ~ {
        Extract_StochScores_DR(
          split_data = .x, 
          fit_nuisance = cfg$fit_nuis, 
          fit_train = cfg$fit_train,
          CV_nusiance = cvs_nuisance, 
          CV_model = cvs_model,
          conformal_type = ctype,    
          cap = FALSE                  
        )
      }, .options = furrr_options(seed = 2026, packages = c("SuperLearner","ranger","BART","BAS","dplyr")))
      

      # bind them all together and tag them with the Repetition ID, Method, and Case.
      method_scores <- bind_rows(res_list, .id = "Rep") %>%
        mutate(
          Rep = as.integer(Rep),
          Method = m,
          Case = case
        ) %>%
        relocate(Case, Method, Rep) # Just formatting the columns nicely
      
      all_scores[[paste0("case", case, "_", m)]] <- method_scores
    }
  }
  bind_rows(all_scores)
}
## ---- make results directories ----
if (!dir.exists("results_stoch")) dir.create("results_stoch", showWarnings = FALSE)

## ---- read Slurm array ID to choose scenario ----
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
cat("SLURM_ARRAY_TASK_ID =", task_id, "\n")

config <- expand.grid(
  case   = 1:36, 
  ctype  = "res", # "clip"
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>% arrange(case, ctype)

if (task_id < 1 || task_id > nrow(config)) stop("task_id out of range.")

this_case  <- config$case[task_id]
this_ctype <- config$ctype[task_id]
scenario <- glue("case{this_case}_{this_ctype}")
cat("Running Stochastic Extraction for scenario:", scenario, "\n")

slurm_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
n_workers <- max(1, slurm_cores - 1)
future::plan(future::multicore, workers = n_workers)
cat("Running on", n_workers, "parallel workers using multicore.\n")

grid_out_scores <- run_grid_stoch(
  cases = this_case,
  ctype = this_ctype,
  cvs_model = TRUE, 
  cvs_nuisance = TRUE,
  methods = c("ID_reg","ID_rf","ID_sl","ID_bart","ID_bma",
              "ED_reg","ED_rf","ED_sl","ED_bart","ED_bma"),
  workers = n_workers
)

rdata_file <- file.path("results_stoch", glue("stoch_scores_{scenario}.RData"))
save(grid_out_scores, file = rdata_file)
