# pre-code
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
RNGkind("L'Ecuyer-CMRG")
library(tidyverse)
library(parallel)
library(ggplot2)
library(gridExtra) # for plot arrange
n_rep <- 1000
n_cores <- 1
seed_base <- 2026
source("0_functions_DGP.R")
#
# Setup ------------------------------------------------------------------------
# sample size
n <- c(50,100,200) # RCT size: 200, 400 
ne <- 800 # ED size
m <- 100 # testing size: 500

p <- 2 # number of covariates

etaj <- c(0.3,0.8) # adjust covariate shift degree (S~X model)
sd_noise_spec <- c("const:0.3","const:0.6","hetero:ITE") #c("const:0.3", "const:0.5", "const:0.6", "const:0.8", "hetero:ITE") #noise in outcome (RCT+EC) #sd_noise <- 0.5

borrow_trt <- TRUE #c(TRUE,FALSE)
linear_structure <- c(TRUE,FALSE)

# hidden bias
prop_unbias <- 0.5
mag_bias <- 0.8 #0.4 # 0 #0.8

setup <- expand_grid(
  borrow_trt  = borrow_trt,
  n = n,
  ne = ne,
  m = m,
  p = p,
  etaj = etaj,
  sd_noise_spec = sd_noise_spec,
  prop_unbias = prop_unbias,
  mag_bias = mag_bias,
  linear_structure = linear_structure
) %>%
  mutate(
    N = n+ne, # pooled sample size in stage 1
    prop_unbias = ifelse(mag_bias == 0, 0, prop_unbias),
    case = row_number(),
    ite_type = ifelse(linear_structure, "linear", "nonlinear")) %>% 
  relocate(case, borrow_trt, ite_type, .before = everything())

if (mag_bias == 0){
  # X-shift only
  out_dir <- str_glue("data")
}else{
  # Add Y-shift
  out_dir <- str_glue("data_bias/mag_bias_{mag_bias}")
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
save(setup, file = file.path(out_dir, "setup.RData"))

# Stage 1 Data -----------------------------------------------------------------
for (i_case in which(setup$sd_noise_spec=="hetero:ITE")) { #seq_len(nrow(setup))
  with(setup[i_case, ], {
    message(sprintf("[case %s] Stage 1 data (RCT+EC)", case))
    ite_fun <- get_ite_fun(ite_type)
    
    # ---------- 1) simdata_main ----------
    set.seed(seed_base)
    
    simdata_main <- mclapply(seq_len(n_rep), function(rep_id) {
      rct <- generate_main(
        n = n,
        main_only = TRUE,
        ite_fun = ite_fun,
        sd_noise_spec = sd_noise_spec,
        seed = seed_base + rep_id
      )
      ed <- generate_main(
        n = ne,
        main_only = FALSE,
        ite_fun = ite_fun,
        
        etaj  = etaj,
        borrow_trt = borrow_trt,
        
        sd_noise_spec = sd_noise_spec,
        seed = seed_base + rep_id + 1
      )
      main_item <- list(
        X = rbind(rct$X, ed$X),
        S = c(rct$S, ed$S),
        A = c(rct$A, ed$A),
        mu0 = c(rct$mu0, ed$mu0),
        mu1 = c(rct$mu1, ed$mu1),
        ITE = c(rct$ITE, ed$ITE),
        Y0 = c(rct$Y0, ed$Y0),
        Y1 = c(rct$Y1, ed$Y1),
        R2_rct = rct$R2,
        R2_ed = ed$R2,
        sd_noise_spec = sd_noise_spec
      )
      main_item
    }, mc.cores = n_cores)
    save(simdata_main, file = file.path(out_dir, str_glue("simdata_main_{case}.RData")))
    # save(simdata_main, file = str_glue("data_bias/simdata_main_{case}.RData"))
    
    # ---------- 2) simdata (observed data, bias/noise options) ----------
    simdata <- mclapply(seq_len(n_rep), function(rep_id) {
      generate_bias_from_main(
        main_item       = simdata_main[[rep_id]],
        borrow_trt      = borrow_trt,
        mag_bias        = mag_bias,
        prop_unbias     = prop_unbias,
        sd_noise_spec = sd_noise_spec,
        seed            = seed_base + rep_id + 2
      )
    }, mc.cores = n_cores)
    save(simdata, file = file.path(out_dir, str_glue("simdata_{case}.RData")))
    # save(simdata, file = str_glue("data_bias/simdata_{case}.RData"))
  })
}

# Testing Data -----------------------------------------------------------------
rename_Y_to_Yt <- function(obj) {
  if (!is.null(obj$Y1)) { obj$Y1t <- obj$Y1; obj$Y1 <- NULL }
  if (!is.null(obj$Y0)) { obj$Y0t <- obj$Y0; obj$Y0 <- NULL }
  if (!is.null(obj$X)) { obj$Xt <- obj$X; obj$X <- NULL }
  obj
}
for (i_case in which(setup$sd_noise_spec=="hetero:ITE")) { # seq_len(nrow(setup))
  with(setup[i_case, ], {
    message(sprintf("[case %s] start generating TEST (S=1-only)", case))
    ite_fun <- get_ite_fun(ite_type)
    simdata_test <- mclapply(seq_len(n_rep), function(rep_id) {
      generate_main(
        n = m,
        main_only = TRUE,
        ite_fun = ite_fun,
        sd_noise_spec = sd_noise_spec,
        seed = seed_base + rep_id + 3
      )
    }, mc.cores = n_cores)
    simdata_test <- lapply(simdata_test, rename_Y_to_Yt)
    save(simdata_test, file = file.path(out_dir, str_glue("simdata_test_{case}.RData")))
    #save(simdata_test, file = str_glue("data_bias/simdata_test_{case}.RData"))
  })
}

