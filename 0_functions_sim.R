subset_by_S <- function(main, full, group = 1) {
  keep <- main$S == group
  list(
    S     = main$S[keep],
    X     = main$X[keep, , drop = FALSE],
    A     = main$A[keep],
    Y0    = main$Y0[keep],
    Y1    = main$Y1[keep],
    mu0   = main$mu0[keep],
    mu1   = main$mu1[keep],
    
    Y     = full$Y[keep],
    Y00   = full$Y00[keep],
    Y11   = full$Y11[keep],
    Ynull = full$Ynull[keep]
  )
}

df_X <- function(data){
  stopifnot(!is.null(data$X))
  X_df <- as.data.frame(data$X)
  colnames(X_df) <- paste0("X", seq_len(ncol(X_df)))
  data$X <- NULL
  cbind(data, X_df)
}

split_data <- function(data, a, ed, use_ed, borrow_trt,
                       split_train = 0.5, split_nuisance = 0.2,
                       method = "DR"){
  dat <- list(
    X  = data$X[data$A == a, ],
    Y  = data$Y[data$A == a],
    Y1 = data$Y1[data$A == a],
    Y0 = data$Y0[data$A == a],
    S = data$S[data$A == a]
  )
  # define id seq for RCT data
  id_rct <- seq_along(dat$Y)
  
  # split the subset into train and calibration (or nuisance)
  if (method == "DR"){
    n_nuisance <- ceiling(split_nuisance * length(dat$Y))
    n_train <- ceiling(split_train * (length(dat$Y) - n_nuisance))
    n_cal <- length(dat$Y) - n_train - n_nuisance
    
    split_id <- split(id_rct,
                      sample(c(rep(1, n_train), rep(2, n_cal), rep(3, n_nuisance))))
    id_train <- split_id[[1]]
    id_cal <- split_id[[2]]
    id_nuisance <- split_id[[3]]
    
    dat_calib <- tibble(Y = dat$Y[id_cal], X = dat$X[id_cal, , drop = FALSE], S = dat$S[id_cal])
    dat_nuisance <- tibble(Y = dat$Y[id_nuisance], X = dat$X[id_nuisance, , drop = FALSE], S = dat$S[id_nuisance])
  } else {
    n_train <- ceiling(split_train * length(dat$Y))
    n_cal <- length(dat$Y) - n_train
    
    split_id <- split(id_rct, sample(c(rep(1, n_train), rep(2, n_cal))))
    id_train <- split_id[[1]]
    id_cal <- split_id[[2]]
    
    dat_calib <- tibble(Y = dat$Y[id_cal], X = dat$X[id_cal, , drop = FALSE], S = dat$S[id_cal])
    dat_nuisance <- NULL
  }
  
  if (!use_ed){
    dat_train <- tibble(Y = dat$Y[id_train], X = dat$X[id_train, , drop = FALSE], S = dat$S[id_train])
  } else {
    if (a == 0){
      ed_train  <- tibble(Y = ed$Y[ed$A == a], X = ed$X[ed$A == a, , drop = FALSE], S = ed$S[ed$A == a])
      rct_train <- tibble(Y = dat$Y[id_train], X = dat$X[id_train, , drop = FALSE], S = dat$S[id_train])
      dat_train <- bind_rows(rct_train, ed_train)
    } else {
      if (isTRUE(borrow_trt)){
        ed_train  <- tibble(Y = ed$Y[ed$A == a], X = ed$X[ed$A == a, , drop = FALSE], S = ed$S[ed$A == a])
        rct_train <- tibble(Y = dat$Y[id_train], X = dat$X[id_train, , drop = FALSE], S = dat$S[id_train])
        dat_train <- bind_rows(rct_train, ed_train)
      } else {
        dat_train <- tibble(Y = dat$Y[id_train], X = dat$X[id_train, , drop = FALSE], S = dat$S[id_train])
      }
    }
  }
  lst(dat_train, dat_calib, dat_nuisance)
}

# non-ties
# rank_pvals <- function(s_cal, s_test){
#   ncal <- length(s_cal)
#   vapply(s_test, function(st) (1 + sum(s_cal >= st)) / (ncal + 1), numeric(1))
# }

rank_pvals <- function(s_cal, s_test, tie=TRUE){
  ncal <- length(s_cal)
  vapply(s_test, function(st) {
    gt <- sum(s_cal < st)
    eq <- sum(s_cal == st)
    u  <- if (tie) runif(1) else 1
    (gt + (1 + eq) * u) / (ncal + 1)
  }, numeric(1))
}

# ConfPval_DR <- function(dat_train, # combined
#                         dat_calib, # combined
#                         dat_nuisance, # combined
#                         X_test, alpha, ps=0.5, 
#                         CV_nusiance = FALSE, # TRUE: use ED+RCT get nuisance para; FALSE: use RCT only to get nuisance para
#                         # CV_ratio = 0.75, # calibration = 0.25
#                         CV_fold = 5,
#                         fit_nuisance=c("reg", "rf", "sl", "bart", "bma"),
#                         fit_train=c("reg", "rf", "sl", "bart", "bma"),
#                         SL_lib_nuisance = SL_LIB_NUIS,
#                         SL_lib_ar = SL_LIB_AR,
#                         CV_model = FALSE,
#                         conformal_type=c("res","clip","CQR"), # only ar can do superlearner
#                         y_star = 0, # for testing points Y_pseudo
#                         return_scores = TRUE,
#                         cap = FALSE # for cf =  <- ax{0,Y_pseu-Y_fitted} or not
# ){
#   # prepare data_all
#   dat_all <- bind_rows(dat_nuisance, dat_train, dat_calib)
#   df_all <- df_X(dat_all) # turn X matrix to be X df for random forest modeling
#   X_cols <- df_all %>% select(starts_with("X")) %>% names()
#   
#   # prepare PS
#   if (!"ps" %in% names(df_all)){
#     if ("S" %in% names(df_all)){
#       df_all$ps <- NA_real_
#       
#       idx_RCT <- which(df_all$S == 1)
#       idx_ED <- which(df_all$S == 0)
#       
#       df_all$ps[idx_RCT] <- ps # defined in argument
#       if (length(idx_ED)>0){
#         df_ED <- df_all[idx_ED, , drop = FALSE]
#         ps_mod <- glm(A ~ .,data = df_ED[, c("A", X_cols), drop = FALSE],family = binomial())
#         df_all$ps[idx_ED] <- predict(ps_mod,newdata = df_ED[, X_cols, drop = FALSE], type = "response")
#       }
#     }else{
#       df_all$ps <- ps
#     }
#   }
#   ps_i <- df_all$ps
#   
#   if (!CV_nusiance){
#     df_nuisance <- df_X(dat_nuisance)
#     # estimate nuisance parameters: ps, hatmu1, hatmu0
#     if (fit_nuisance == "rf"){
#       train1 <- df_nuisance %>%
#         dplyr::filter(A == 1) %>%
#         dplyr::select(Y, dplyr::all_of(X_cols))
#       train0 <- df_nuisance %>%
#         dplyr::filter(A == 0) %>%
#         dplyr::select(Y, dplyr::all_of(X_cols))
#       fit1 <- ranger::ranger(
#         Y ~ ., 
#         data = train1,
#         num.trees = 100,
#         max.depth = 5,
#         seed = 0
#       )
#       fit0 <- ranger::ranger(
#         Y ~ ., 
#         data = train0,
#         num.trees = 100,
#         max.depth = 5,
#         seed = 0
#       )
#       X_all <- df_all[, X_cols, drop = FALSE]
#       hatmu1 <- predict(fit1, data = X_all)$predictions
#       hatmu0 <- predict(fit0, data = X_all)$predictions
#     }else if (fit_nuisance == "reg"){
#       hatmu1 <- lm(Y ~ ., data = df_nuisance %>% filter(A==1) %>% select(Y, all_of(X_cols))) %>% predict(newdata = df_all[, X_cols, drop = FALSE])
#       hatmu0 <- lm(Y ~ ., data = df_nuisance %>% filter(A==0) %>% select(Y, all_of(X_cols))) %>% predict(newdata = df_all[, X_cols, drop = FALSE])
#     }else{
#       df1 <- df_nuisance %>% dplyr::filter(A==1)
#       df0 <- df_nuisance %>% dplyr::filter(A==0)
#       fit1 <- SuperLearner::SuperLearner(
#         Y = df1$Y, X = df1[, X_cols, drop=FALSE],
#         SL.library = SL_lib_nuisance, family = stats::gaussian()
#       )
#       fit0 <- SuperLearner::SuperLearner(
#         Y = df0$Y, X = df0[, X_cols, drop=FALSE],
#         SL.library = SL_lib_nuisance, family = stats::gaussian()
#       )
#       hatmu1 <- as.numeric(predict(fit1, newdata = df_all[, X_cols, drop=FALSE])$pred)
#       hatmu0 <- as.numeric(predict(fit0, newdata = df_all[, X_cols, drop=FALSE])$pred)
#     }
#     # generate pseudo outcomes using DR
#     A <- df_all$A
#     Y <- df_all$Y
#     # df_all$Y_pseudo <- (A-ps)/(ps*(1-ps))*(Y-(hatmu0*(1-A)+hatmu1*A)) + hatmu1-hatmu0
#     df_all$Y_pseudo <- (A-ps_i)/(ps_i*(1-ps_i))*(Y-(hatmu0*(1-A)+hatmu1*A)) + hatmu1-hatmu0
#     
#     # train model on pseudo outcome & get conformal scores on calibration set
#     df_train <- df_all %>% filter(type=="train")
#     df_calib <- df_all %>% filter(type=="calibration")
#   }else{
#     
#     df_model <- df_all %>% dplyr::filter(type %in% c("train","nuisance"))
#     df_calib <- df_all %>% dplyr::filter(type == "calibration")
#     # set.seed(seed)
#     # if (!is.null(CV_ratio) && CV_ratio <= 1) {
#     #   n_model <- nrow(df_model)
#     #   idx_model <- sample(seq_len(n_model), size = floor(CV_ratio * n_model))
#     #   df_model  <- df_model[idx_model, , drop = FALSE]
#     # }
#     
#     X_cols <- df_model %>% dplyr::select(dplyr::starts_with("X")) %>% names()
#     K <- CV_fold
#     folds <- caret::createFolds(df_model$Y, k = K)
#     Y_pseudo_model <- rep(NA_real_, nrow(df_model))
#     Y_pseudo_calib_mat <- matrix(NA_real_, nrow = nrow(df_calib), ncol = K)
#     
#     for (k in seq_len(K)) {
#       fold_idx  <- folds[[k]]
#       train_idx <- setdiff(seq_len(nrow(df_model)), fold_idx)
#       df_train_k <- df_model[train_idx, , drop = FALSE]
#       
#       if (fit_nuisance == "rf") {
#         train1_k <- df_train_k %>%
#           dplyr::filter(A == 1) %>%
#           dplyr::select(Y, dplyr::all_of(X_cols))
#         
#         train0_k <- df_train_k %>%
#           dplyr::filter(A == 0) %>%
#           dplyr::select(Y, dplyr::all_of(X_cols))
#         
#         fit1 <- ranger::ranger(
#           Y ~ .,
#           data      = train1_k,
#           num.trees = 100,
#           max.depth = 5,
#           seed      = 0
#         )
#         fit0 <- ranger::ranger(
#           Y ~ .,
#           data      = train0_k,
#           num.trees = 100,
#           max.depth = 5,
#           seed      = 0
#         )
#         X_model_fold <- df_model[fold_idx, X_cols, drop = FALSE]
#         X_calib_all <- df_calib[, X_cols, drop = FALSE]
#         
#         hatmu1_k <- predict(fit1, data = X_model_fold)$predictions
#         hatmu0_k <- predict(fit0, data = X_model_fold)$predictions
#         hatmu1_calib <- predict(fit1, data = X_calib_all)$predictions
#         hatmu0_calib <- predict(fit0, data = X_calib_all)$predictions
#         
#       } else if (fit_nuisance == "reg") {
#         fit1 <- stats::lm(Y ~ ., data = df_train_k %>%
#                             dplyr::filter(A == 1) %>%
#                             dplyr::select(Y, dplyr::all_of(X_cols)))
#         fit0 <- stats::lm(Y ~ ., data = df_train_k %>%
#                             dplyr::filter(A == 0) %>%
#                             dplyr::select(Y, dplyr::all_of(X_cols)))
#         
#         hatmu1_k <- predict(fit1, newdata = df_model[fold_idx, X_cols, drop = FALSE])
#         hatmu0_k <- predict(fit0, newdata = df_model[fold_idx, X_cols, drop = FALSE])
#         hatmu1_calib <- predict(fit1, newdata = df_calib[, X_cols, drop = FALSE])
#         hatmu0_calib <- predict(fit0, newdata = df_calib[, X_cols, drop = FALSE])
#         
#       } else {  # "sl"
#         df1_k <- dplyr::filter(df_train_k, A == 1)
#         df0_k <- dplyr::filter(df_train_k, A == 0)
#         
#         sl1 <- SuperLearner::SuperLearner(
#           Y = df1_k$Y, X = df1_k[, X_cols, drop = FALSE],
#           SL.library = SL_lib_nuisance, family = stats::gaussian()
#         )
#         sl0 <- SuperLearner::SuperLearner(
#           Y = df0_k$Y, X = df0_k[, X_cols, drop = FALSE],
#           SL.library = SL_lib_nuisance, family = stats::gaussian()
#         )
#         
#         hatmu1_k     <- as.numeric(predict(sl1, newdata = df_model[fold_idx, X_cols, drop = FALSE])$pred)
#         hatmu0_k     <- as.numeric(predict(sl0, newdata = df_model[fold_idx, X_cols, drop = FALSE])$pred)
#         hatmu1_calib <- as.numeric(predict(sl1, newdata = df_calib[, X_cols, drop = FALSE])$pred)
#         hatmu0_calib <- as.numeric(predict(sl0, newdata = df_calib[, X_cols, drop = FALSE])$pred)
#       }
#       A_k <- df_model$A[fold_idx]
#       Y_k <- df_model$Y[fold_idx]
#       ps_k <- df_model$ps[fold_idx]
#       #Y_pseudo_model[fold_idx] <- (A_k - ps) / (ps * (1 - ps)) * (Y_k - (hatmu0_k * (1 - A_k) + hatmu1_k * A_k)) + (hatmu1_k - hatmu0_k)
#       Y_pseudo_model[fold_idx] <- (A_k - ps_k) / (ps_k * (1 - ps_k)) * (Y_k - (hatmu0_k * (1 - A_k) + hatmu1_k * A_k)) + (hatmu1_k - hatmu0_k)
#       
#       A_c <- df_calib$A
#       Y_c <- df_calib$Y
#       ps_c <- df_calib$ps
#       Y_pseudo_calib_mat[, k] <- (A_c - ps_c) / (ps_c * (1 - ps_c)) * (Y_c - (hatmu0_calib * (1 - A_c) + hatmu1_calib * A_c)) + (hatmu1_calib - hatmu0_calib)
#     }
#     df_train <- df_model
#     df_train$Y_pseudo <- Y_pseudo_model
#     df_calib$Y_pseudo <- rowMeans(Y_pseudo_calib_mat)
#   }
#   
#   # prepare functions to get one-sided conformal p-value
#   X_test_df <- as.data.frame(X_test)
#   colnames(X_test_df) <- paste0("X", 1:ncol(X_test_df))
#   y_calib   <- df_calib$Y_pseudo
#   
#   # get conformal scores
#   if (conformal_type %in% c("res","clip")){
#     if(!CV_model){
#       if (fit_train == "reg"){
#         hatmu_pseudo <- lm(Y_pseudo ~ ., data = df_train %>% select(Y_pseudo, all_of(X_cols)))
#         mu_cal <- predict(hatmu_pseudo, newdata = df_calib) # fitted outcome for calibration
#         mu_tst <- predict(hatmu_pseudo, newdata = X_test_df) # fitted outcome for testing
#       } else if (fit_train == "rf"){
#         train_rf <- df_train %>% dplyr::select(Y_pseudo, dplyr::all_of(X_cols))
#         hatmu_pseudo <- ranger::ranger(
#           Y_pseudo ~ ., 
#           data = train_rf,
#           num.trees = 100,
#           max.depth = 5,
#           seed = 0)
#         mu_cal <- predict(hatmu_pseudo, data = df_calib %>% dplyr::select(dplyr::all_of(X_cols)))$predictions # fitted outcome for calibration
#         mu_tst <- predict(hatmu_pseudo, data = X_test_df)$predictions # fitted outcome for testing
#       } else {
#         sl_fit <- SuperLearner::SuperLearner(
#           Y = df_train$Y_pseudo, X = df_train[, X_cols, drop=FALSE],
#           SL.library = SL_lib_ar, family = stats::gaussian()
#         )
#         mu_cal <- as.numeric(predict(sl_fit, newdata = df_calib[, X_cols, drop=FALSE])$pred)
#         mu_tst <- as.numeric(predict(sl_fit, newdata = X_test_df[, X_cols, drop=FALSE])$pred)
#       }
#     }else{ # CV in train
#       K_pred <- CV_fold
#       folds_model <- caret::createFolds(df_train$Y_pseudo, k = K_pred)
#       n_cal <- nrow(df_calib)
#       n_tst <- nrow(X_test_df)
#       mu_cal_mat <- matrix(NA_real_, nrow = n_cal, ncol = K_pred)
#       mu_tst_mat <- matrix(NA_real_, nrow = n_tst, ncol = K_pred)
#       
#       for (k in seq_len(K_pred)) {
#         fold_idx_pred  <- folds_model[[k]]
#         train_idx_pred <- setdiff(seq_len(nrow(df_train)), fold_idx_pred)
#         df_train_k_pred <- df_train[train_idx_pred, , drop = FALSE]
#         
#         if (fit_train == "reg") {
#           fit_k <- lm(Y_pseudo ~ ., data = df_train_k_pred %>% dplyr::select(Y_pseudo, dplyr::all_of(X_cols)))
#           mu_cal_mat[, k] <- predict(fit_k, newdata = df_calib)
#           mu_tst_mat[, k] <- predict(fit_k, newdata = X_test_df)
#           
#         } else if (fit_train == "rf") {
#           train_rf_k <- df_train_k_pred %>% dplyr::select(Y_pseudo, dplyr::all_of(X_cols))
#           fit_k <- ranger::ranger(Y_pseudo ~ ., data = train_rf_k,
#             num.trees = 100, max.depth = 5, seed = 0)
#           mu_cal_mat[, k] <- predict(fit_k, data = df_calib %>% dplyr::select(dplyr::all_of(X_cols)))$predictions
#           mu_tst_mat[, k] <- predict(fit_k, data = X_test_df)$predictions
#           
#         } else { # "sl"
#           sl_fit_k <- SuperLearner::SuperLearner(
#             Y = df_train_k_pred$Y_pseudo,
#             X = df_train_k_pred[, X_cols, drop = FALSE],
#             SL.library = SL_lib_ar,
#             family = stats::gaussian()
#           )
#           mu_cal_mat[, k] <- as.numeric(predict(sl_fit_k, newdata = df_calib[, X_cols, drop = FALSE])$pred)
#           mu_tst_mat[, k] <- as.numeric(predict(sl_fit_k, newdata = X_test_df[, X_cols, drop = FALSE])$pred)
#         }
#         
#       }
#       mu_cal <- rowMeans(mu_cal_mat)
#       mu_tst <- rowMeans(mu_tst_mat)
#     }
#     
#     if (conformal_type == "res"){
#       if(cap){
#         s_cal <- pmax(0, y_calib - mu_cal)
#         s_tst <- pmax(0, y_star - mu_tst)
#       }else{
#         s_cal <- (y_calib - mu_cal)
#         s_tst <- (y_star - mu_tst)
#       } 
#     }else if (conformal_type == "clip"){
#       M <- 100
#       if(cap){
#         s_cal <- pmax(0, M * (y_calib > 0) - mu_cal)
#         s_tst <- pmax(0, M * (y_star > 0) - mu_tst)
#       }else{
#         s_cal <- M * (y_calib > 0) - mu_cal
#         s_tst <- M * (y_star > 0) - mu_tst
#       }   
#     }
#     pvals <- rank_pvals(s_cal, s_tst) 
#     # q_calib <- predict(hatmu_pseudo, newdata = df_calib)
#     # cs_calib <- abs(df_calib$Y_pseudo - q_calib)
#   }
#   
#   # if (conformal_type == "CQR"){
#   #   tau_side <- if (side == "upper") 1 - alpha else alpha
#   #   if (fit_train == "reg"){
#   #     qr <- quantreg::rq(Y_pseudo ~ ., tau = tau_side, data = df_train %>% select(Y_pseudo, all_of(X_cols)))
#   #     q_cal <- as.numeric(predict(qr, newdata = df_calib))
#   #     q_tst <- as.numeric(predict(qr, newdata = X_test_df))
#   #   }
#   #   if (fit_train == "rf"){
#   #     qr <- grf::quantile_forest(
#   #       df_train[, X_cols],
#   #       df_train$Y_pseudo,
#   #       quantiles = tau_side
#   #     )
#   #     q_cal <- as.numeric(predict(qr, df_calib[, X_cols])$predictions)
#   #     q_tst <- as.numeric(predict(qr, X_test_df)$predictions)
#   #   }
#   #   if (fit_train == "sl"){
#   #     qSL <- .fit_quantile_dSL(df_train, X_cols, tau = tau_side, V = 3)
#   #     
#   #     q_cal <- qSL$predict(df_calib[, X_cols, drop = FALSE])
#   #     q_tst <- qSL$predict(X_test_df[, X_cols, drop = FALSE])
#   #   }
#   #   if (side == "upper") {
#   #     if(cap){
#   #       s_cal <- pmax(0, y_calib - q_cal)
#   #       s_tst <- pmax(0, y_star - q_tst)
#   #     }else{
#   #       s_cal <- (y_calib - q_cal)
#   #       s_tst <- (y_star - q_tst)
#   #     }
#   #   } else { # "lower"
#   #     if(cap){
#   #       s_cal <- pmax(0, q_cal - y_calib)
#   #       s_tst <- pmax(0, q_tst - y_star)
#   #     }else{
#   #       s_cal <- (q_cal - y_calib)
#   #       s_tst <- (q_tst - y_star)
#   #     }
#   #   }
#   #   pvals <- rank_pvals(s_cal, s_tst) 
#   # }
#   
#   
#   # if we want to add cs values in output
#   p_raw <- pvals
#   p_BH <- p.adjust(p_raw, method = "BH")
#   p_Bon <- p.adjust(p_raw, method = "bonferroni")
#   out <- tibble::tibble(
#     p_value = p_raw,
#     p_BH = p_BH,    # BH-adjusted
#     p_Bon = p_Bon    # Bonferroni-adjusted
#   )
#   
#   if (return_scores) {
#     # for calibration set
#     out$s_calib <- list(s_cal) # conformal score
#     out$ITE_pseudo_calib <- list(df_calib$Y_pseudo) # pseudo ITE
#     out$ITE_fitted_calib <- list(mu_cal)
#     
#     # for testing set
#     out$s_test  <- s_tst
#     out$ITE_fitted_test <- as.numeric(mu_tst)
#     
#     out$alpha <- alpha
#   }
#   out
# }
fit_predict_gaussian <- function(df_train, df_pred, y_col, x_cols,
                                 method = c("reg","rf","sl","bart","bma"),
                                 SL_lib = NULL,
                                 rf_num_trees = 100,
                                 rf_max_depth = 5,
                                 seed = 0) {
  method <- match.arg(method)
  Xtr <- df_train[, x_cols, drop = FALSE]
  Xpr <- df_pred[,  x_cols, drop = FALSE]
  y   <- df_train[[y_col]]
  
  if (method == "reg") {
    fit <- stats::lm(stats::as.formula(paste0(y_col, " ~ .")),
                     data = df_train[, c(y_col, x_cols), drop = FALSE])
    return(as.numeric(stats::predict(fit, newdata = df_pred[, x_cols, drop = FALSE])))
    
  } else if (method == "rf") {
    fit <- ranger::ranger(
      formula = stats::as.formula(paste0(y_col, " ~ .")),
      data    = df_train[, c(y_col, x_cols), drop = FALSE],
      num.trees = rf_num_trees,
      max.depth = rf_max_depth,
      seed = seed
    )
    return(as.numeric(predict(fit, data = Xpr)$predictions))
    
  } else if (method == "sl") {
    if (is.null(SL_lib)) stop("SL_lib must be provided when method = 'sl'.")
    fit <- SuperLearner::SuperLearner(
      Y = y,
      X = Xtr,
      SL.library = SL_lib,
      family = stats::gaussian()
    )
    return(as.numeric(predict(fit, newdata = Xpr)$pred))
    
  } else if (method == "bart") {
    # BART wants matrices
    Xtr_m <- as.matrix(Xtr)
    Xpr_m <- as.matrix(Xpr)
    fit <- BART::gbart(
      x.train = Xtr_m,
      y.train = y,
      x.test  = Xpr_m,
      ndpost = 200L,
      ntree  = 50L
    )
    return(as.numeric(colMeans(fit$yhat.test)))
    
  } else if (method == "bma") {
    # Bayesian model averaging for linear regression (continuous outcomes)
    # BAS::bas.lm supports formula interface; prediction is easy.
    if (!requireNamespace("BAS", quietly = TRUE)) {
      stop("Package 'BAS' is required for method='bma'. Install with install.packages('BAS').")
    }
    df_fit <- df_train[, c(y_col, x_cols), drop = FALSE]
    fit <- BAS::bas.lm(
      formula = stats::as.formula(paste0(y_col, " ~ .")),
      data    = df_fit,
      method  = "BAS",     # robust default sampler
      prior   = "g-prior", # common default
      modelprior = BAS::uniform()  # uniform over models
    )
    pr <- predict(fit, newdata = df_pred[, x_cols, drop = FALSE], estimator = "BMA")
    return(as.numeric(pr$fit))
  }
}

ConfPval_DR <- function(dat_train, # combined
                        dat_calib, # combined
                        dat_nuisance, # combined
                        X_test, alpha, ps=0.5, 
                        CV_nusiance = FALSE, # TRUE: use ED+RCT get nuisance para; FALSE: use RCT only to get nuisance para
                        # CV_ratio = 0.75, # calibration = 0.25
                        CV_fold = 5,
                        fit_nuisance=c("reg", "rf", "sl", "bart", "bma"),
                        fit_train=c("reg", "rf", "sl", "bart", "bma"),
                        SL_lib_nuisance = SL_LIB_NUIS,
                        SL_lib_ar = SL_LIB_AR,
                        CV_model = FALSE,
                        conformal_type=c("res","clip","CQR"), # only ar can do superlearner
                        y_star = 0, # for testing points Y_pseudo
                        return_scores = TRUE,
                        cap = FALSE # for cf =  <- ax{0,Y_pseu-Y_fitted} or not
){
  # prepare data_all
  dat_all <- bind_rows(dat_nuisance, dat_train, dat_calib)
  df_all <- df_X(dat_all) # turn X matrix to be X df for random forest modeling
  X_cols <- df_all %>% select(starts_with("X")) %>% names()
  
  # prepare PS
  if (!"ps" %in% names(df_all)){
    if ("S" %in% names(df_all)){
      df_all$ps <- NA_real_
      
      idx_RCT <- which(df_all$S == 1)
      idx_ED <- which(df_all$S == 0)
      
      df_all$ps[idx_RCT] <- ps # defined in argument
      if (length(idx_ED)>0){
        df_ED <- df_all[idx_ED, , drop = FALSE]
        ps_mod <- glm(A ~ .,data = df_ED[, c("A", X_cols), drop = FALSE],family = binomial())
        df_all$ps[idx_ED] <- predict(ps_mod,newdata = df_ED[, X_cols, drop = FALSE], type = "response")
      }
    }else{
      df_all$ps <- ps
    }
  }
  ps_i <- df_all$ps
  
  if (!CV_nusiance){
    df_nuisance <- df_X(dat_nuisance)
    # estimate nuisance parameters: ps, hatmu1, hatmu0
    train1 <- df_nuisance %>%
      dplyr::filter(A == 1) %>%
      dplyr::select(Y, dplyr::all_of(X_cols))
    train0 <- df_nuisance %>%
      dplyr::filter(A == 0) %>%
      dplyr::select(Y, dplyr::all_of(X_cols))
   
     hatmu1 <- fit_predict_gaussian(
      df_train = train1,
      df_pred  = df_all,
      y_col    = "Y",
      x_cols   = X_cols,
      method   = fit_nuisance,
      SL_lib   = SL_lib_nuisance
    )
    hatmu0 <- fit_predict_gaussian(
      df_train = train0,
      df_pred  = df_all,
      y_col    = "Y",
      x_cols   = X_cols,
      method   = fit_nuisance,
      SL_lib   = SL_lib_nuisance
    )
    
   
    # generate pseudo outcomes using DR
    A <- df_all$A
    Y <- df_all$Y
    # df_all$Y_pseudo <- (A-ps)/(ps*(1-ps))*(Y-(hatmu0*(1-A)+hatmu1*A)) + hatmu1-hatmu0
    df_all$Y_pseudo <- (A-ps_i)/(ps_i*(1-ps_i))*(Y-(hatmu0*(1-A)+hatmu1*A)) + hatmu1-hatmu0
    
    # train model on pseudo outcome & get conformal scores on calibration set
    df_train <- df_all %>% filter(type=="train")
    df_calib <- df_all %>% filter(type=="calibration")
  }else{
    
    df_model <- df_all %>% dplyr::filter(type %in% c("train","nuisance"))
    df_calib <- df_all %>% dplyr::filter(type == "calibration")
    # set.seed(seed)
    # if (!is.null(CV_ratio) && CV_ratio <= 1) {
    #   n_model <- nrow(df_model)
    #   idx_model <- sample(seq_len(n_model), size = floor(CV_ratio * n_model))
    #   df_model  <- df_model[idx_model, , drop = FALSE]
    # }
    
    X_cols <- df_model %>% dplyr::select(dplyr::starts_with("X")) %>% names()
    K <- CV_fold
    folds <- caret::createFolds(df_model$Y, k = K)
    Y_pseudo_model <- rep(NA_real_, nrow(df_model))
    Y_pseudo_calib_mat <- matrix(NA_real_, nrow = nrow(df_calib), ncol = K)
    
    for (k in seq_len(K)) {
      fold_idx  <- folds[[k]]
      train_idx <- setdiff(seq_len(nrow(df_model)), fold_idx)
      df_train_k <- df_model[train_idx, , drop = FALSE]
      
      train1_k <- df_train_k %>%
        dplyr::filter(A == 1) %>%
        dplyr::select(Y, dplyr::all_of(X_cols))
      
      train0_k <- df_train_k %>%
        dplyr::filter(A == 0) %>%
        dplyr::select(Y, dplyr::all_of(X_cols))
      
      # predictions for held-out fold rows
      df_fold_pred <- df_model[fold_idx, , drop = FALSE]
      hatmu1_k <- fit_predict_gaussian(
        df_train = train1_k, df_pred = df_fold_pred,
        y_col = "Y", x_cols = X_cols,
        method = fit_nuisance, SL_lib = SL_lib_nuisance
      )
      hatmu0_k <- fit_predict_gaussian(
        df_train = train0_k, df_pred = df_fold_pred,
        y_col = "Y", x_cols = X_cols,
        method = fit_nuisance, SL_lib = SL_lib_nuisance
      )
      # predictions for calibration rows
      hatmu1_calib <- fit_predict_gaussian(
        df_train = train1_k, df_pred = df_calib,
        y_col = "Y", x_cols = X_cols,
        method = fit_nuisance, SL_lib = SL_lib_nuisance
      )
      hatmu0_calib <- fit_predict_gaussian(
        df_train = train0_k, df_pred = df_calib,
        y_col = "Y", x_cols = X_cols,
        method = fit_nuisance, SL_lib = SL_lib_nuisance
      )
      
      A_k <- df_model$A[fold_idx]
      Y_k <- df_model$Y[fold_idx]
      ps_k <- df_model$ps[fold_idx]
      #Y_pseudo_model[fold_idx] <- (A_k - ps) / (ps * (1 - ps)) * (Y_k - (hatmu0_k * (1 - A_k) + hatmu1_k * A_k)) + (hatmu1_k - hatmu0_k)
      Y_pseudo_model[fold_idx] <- (A_k - ps_k) / (ps_k * (1 - ps_k)) * (Y_k - (hatmu0_k * (1 - A_k) + hatmu1_k * A_k)) + (hatmu1_k - hatmu0_k)
      
      A_c <- df_calib$A
      Y_c <- df_calib$Y
      ps_c <- df_calib$ps
      Y_pseudo_calib_mat[, k] <- (A_c - ps_c) / (ps_c * (1 - ps_c)) * (Y_c - (hatmu0_calib * (1 - A_c) + hatmu1_calib * A_c)) + (hatmu1_calib - hatmu0_calib)
    }
    df_train <- df_model
    df_train$Y_pseudo <- Y_pseudo_model
    df_calib$Y_pseudo <- rowMeans(Y_pseudo_calib_mat)
  }
  
  # prepare functions to get one-sided conformal p-value
  X_test_df <- as.data.frame(X_test)
  colnames(X_test_df) <- paste0("X", 1:ncol(X_test_df))
  y_calib   <- df_calib$Y_pseudo
  
  # get conformal scores
  if (conformal_type %in% c("res","clip")){
    if(!CV_model){
      mu_cal <- fit_predict_gaussian(
        df_train = df_train,
        df_pred  = df_calib,
        y_col    = "Y_pseudo",
        x_cols   = X_cols,
        method   = fit_train,
        SL_lib   = SL_lib_ar
      )
      mu_tst <- fit_predict_gaussian(
        df_train = df_train,
        df_pred  = X_test_df,
        y_col    = "Y_pseudo",
        x_cols   = X_cols,
        method   = fit_train,
        SL_lib   = SL_lib_ar
      )
    }else{ # CV in train
      K_pred <- CV_fold
      folds_model <- caret::createFolds(df_train$Y_pseudo, k = K_pred)
      n_cal <- nrow(df_calib)
      n_tst <- nrow(X_test_df)
      mu_cal_mat <- matrix(NA_real_, nrow = n_cal, ncol = K_pred)
      mu_tst_mat <- matrix(NA_real_, nrow = n_tst, ncol = K_pred)
      
      for (k in seq_len(K_pred)) {
        fold_idx_pred  <- folds_model[[k]]
        train_idx_pred <- setdiff(seq_len(nrow(df_train)), fold_idx_pred)
        df_train_k_pred <- df_train[train_idx_pred, , drop = FALSE]
        
        mu_cal_mat[, k] <- fit_predict_gaussian(
          df_train = df_train_k_pred,
          df_pred  = df_calib,
          y_col    = "Y_pseudo",
          x_cols   = X_cols,
          method   = fit_train,
          SL_lib   = SL_lib_ar
        )
        
        mu_tst_mat[, k] <- fit_predict_gaussian(
          df_train = df_train_k_pred,
          df_pred  = X_test_df,
          y_col    = "Y_pseudo",
          x_cols   = X_cols,
          method   = fit_train,
          SL_lib   = SL_lib_ar
        )
        
      }
      mu_cal <- rowMeans(mu_cal_mat)
      mu_tst <- rowMeans(mu_tst_mat)
    }
    
    if (conformal_type == "res"){
      if(cap){
        s_cal <- pmax(0, y_calib - mu_cal)
        s_tst <- pmax(0, y_star - mu_tst)
      }else{
        s_cal <- (y_calib - mu_cal)
        s_tst <- (y_star - mu_tst)
      } 
    }else if (conformal_type == "clip"){
      M <- 100
      if(cap){
        s_cal <- pmax(0, M * (y_calib > 0) - mu_cal)
        s_tst <- pmax(0, M * (y_star > 0) - mu_tst)
      }else{
        s_cal <- M * (y_calib > 0) - mu_cal
        s_tst <- M * (y_star > 0) - mu_tst
      }   
    }
    pvals <- rank_pvals(s_cal, s_tst) 
    # q_calib <- predict(hatmu_pseudo, newdata = df_calib)
    # cs_calib <- abs(df_calib$Y_pseudo - q_calib)
  }
  
  # if (conformal_type == "CQR"){
  #   tau_side <- if (side == "upper") 1 - alpha else alpha
  #   if (fit_train == "reg"){
  #     qr <- quantreg::rq(Y_pseudo ~ ., tau = tau_side, data = df_train %>% select(Y_pseudo, all_of(X_cols)))
  #     q_cal <- as.numeric(predict(qr, newdata = df_calib))
  #     q_tst <- as.numeric(predict(qr, newdata = X_test_df))
  #   }
  #   if (fit_train == "rf"){
  #     qr <- grf::quantile_forest(
  #       df_train[, X_cols],
  #       df_train$Y_pseudo,
  #       quantiles = tau_side
  #     )
  #     q_cal <- as.numeric(predict(qr, df_calib[, X_cols])$predictions)
  #     q_tst <- as.numeric(predict(qr, X_test_df)$predictions)
  #   }
  #   if (fit_train == "sl"){
  #     qSL <- .fit_quantile_dSL(df_train, X_cols, tau = tau_side, V = 3)
  #     
  #     q_cal <- qSL$predict(df_calib[, X_cols, drop = FALSE])
  #     q_tst <- qSL$predict(X_test_df[, X_cols, drop = FALSE])
  #   }
  #   if (side == "upper") {
  #     if(cap){
  #       s_cal <- pmax(0, y_calib - q_cal)
  #       s_tst <- pmax(0, y_star - q_tst)
  #     }else{
  #       s_cal <- (y_calib - q_cal)
  #       s_tst <- (y_star - q_tst)
  #     }
  #   } else { # "lower"
  #     if(cap){
  #       s_cal <- pmax(0, q_cal - y_calib)
  #       s_tst <- pmax(0, q_tst - y_star)
  #     }else{
  #       s_cal <- (q_cal - y_calib)
  #       s_tst <- (q_tst - y_star)
  #     }
  #   }
  #   pvals <- rank_pvals(s_cal, s_tst) 
  # }
  
  
  # if we want to add cs values in output
  p_raw <- pvals
  p_BH <- p.adjust(p_raw, method = "BH")
  p_Bon <- p.adjust(p_raw, method = "bonferroni")
  out <- tibble::tibble(
    p_value = p_raw,
    p_BH = p_BH,    # BH-adjusted
    p_Bon = p_Bon    # Bonferroni-adjusted
  )
  
  if (return_scores) {
    # for calibration set
    out$s_calib <- list(s_cal) # conformal score
    out$ITE_pseudo_calib <- list(df_calib$Y_pseudo) # pseudo ITE
    out$ITE_fitted_calib <- list(mu_cal)
    
    # for testing set
    out$s_test  <- s_tst
    out$ITE_fitted_test <- as.numeric(mu_tst)
    
    out$alpha <- alpha
  }
  out
}

get_Pval_DR <- function(split_data, test_data, alpha = alpha_overall, ps = 0.5,
                        fit_nuisance = "rf", fit_train = "rf", conformal_type = "clip", CV_nusiance, CV_model,
                        SL_lib_ar = SL_LIB_AR,
                        SL_lib_nuisance = SL_LIB_NUIS,
                        cap = FALSE){
  dat_nuisance1 <- split_data$split_dat1$dat_nuisance %>% mutate(A = 1)
  dat_nuisance0 <- split_data$split_dat0$dat_nuisance %>% mutate(A = 0)
  dat_train1 <- split_data$split_dat1$dat_train %>% mutate(A = 1)
  dat_train0 <- split_data$split_dat0$dat_train %>% mutate(A = 0)
  dat_calib1 <- split_data$split_dat1$dat_calib %>% mutate(A = 1)
  dat_calib0 <- split_data$split_dat0$dat_calib %>% mutate(A = 0)
  
  dat_nuisance <- bind_rows(dat_nuisance1, dat_nuisance0) %>% mutate(type = "nuisance")
  dat_calib <- bind_rows(dat_calib1, dat_calib0) %>% mutate(type = "calibration")
  dat_train <- bind_rows(dat_train1, dat_train0) %>% mutate(type = "train")
  
  X_test <- test_data$Xt
  
  ConfPval_DR(
    dat_train = dat_train, 
    dat_calib = dat_calib, 
    dat_nuisance = dat_nuisance,
    X_test = X_test,
    alpha = alpha,
    ps = ps,
    fit_nuisance = fit_nuisance,
    fit_train = fit_train,
    SL_lib_nuisance = SL_lib_nuisance,
    SL_lib_ar = SL_lib_ar,
    CV_nusiance = CV_nusiance,
    CV_model = CV_model,
    conformal_type = conformal_type,
    cap = cap)
  
}

auc_from_pvals <- function(pval, ite_true) {
  y <- ite_true > 0  # beneficial
  m <- sum(y)
  n <- sum(!y)
  if (m == 0 || n == 0) return(NA_real_)   # undefined if only one class
  
  s <- -pval
  r <- rank(s, ties.method = "average")
  (sum(r[y]) - m*(m + 1)/2) / (m * n)
}

evaluate_metrics_pval <- function(Pval_ite, test_ite_true, gamma,
                                  alpha_FDR = 0.05,
                                  FDR_method = c("None","BH","Bonferroni","mLORD"),
                                  # mlord paras:
                                  w0 = 0.14, b0 = 0.01, decay = 0.3,
                                  xi_mod = NULL){
  # handle NA
  ok <- is.finite(Pval_ite) & is.finite(test_ite_true)
  p  <- Pval_ite[ok]
  t  <- test_ite_true[ok]
  
  H0 <- t <= 0
  H1 <- !H0
  
  if (FDR_method == "None"){
    accept <- p > gamma
    reject <- p <= gamma
    rule <- "No FDR control"
  } else if (FDR_method == "BH") {
    qvals  <- p.adjust(p, method = "BH")
    reject <- (qvals <= alpha_FDR)
    accept <- !reject
    rule <- sprintf("BH (alpha = %.3f)", alpha_FDR)
  } else if (FDR_method == "Bonferroni"){
    m <- length(p)
    qvals  <- p.adjust(p, method = "bonferroni")
    reject <- (qvals <= alpha_FDR) 
    accept <- !reject
    rule   <- sprintf("Bonferroni (alpha = %.3f, m = %d)", alpha_FDR, length(p))
  } else if (FDR_method == "mLORD"){
    if (is.null(xi_mod)) {
      xi_mod <- generate_xi_mod(length(p), alpha = alpha_FDR, b0 = b0, decay = decay)
    }
    res_mlord <- run_modified_lord(pvals = p, alpha = alpha_FDR, w0 = w0, b0 = b0, decay = decay, xi = xi_mod)
    reject <- as.logical(res_mlord$R_i)
    accept <- !reject
    rule <- sprintf("mLORD (alpha = %.3f, w0 = %.3g, b0 = %.3g, decay = %.2f)",
                    alpha_FDR, w0, b0, decay)
  }
  
  R <- sum(accept) # accept (accept H0)
  D <- sum(reject) # discover (reject H0)
  N <- sum(H0)
  P <- sum(H1)
  
  FP <- sum(reject & H0)
  TP <- sum(reject & H1)
  TN <- sum(accept & H0)
  FN <- sum(accept & H1)
  
  FDP_val <- if ( D == 0 ) 0 else FP/D # how many discoveries are beneficial
  power_val <- if (P == 0) NA_real_ else TP/P # given nonbeneficial and being removed
  
  tibble::tibble(
    sensitivity = if (P == 0) NA_real_ else TP/P, # (reject H0 & bene) / bene 
    specificity = if (N == 0) NA_real_ else TN/N, # (accept H0 & nonbene) / nonbene
    precision = if ((TP+FP) == 0) 0 else TP/(TP+FP), # PPV mean(bene[accept]),
    accuracy = mean((accept & H0) | (reject & H1)), 
    prob_nonbene_select = if (N == 0) NA_real_ else FP/N, # prob of accept / nonbene
    auc = auc_from_pvals(p, t),
    
    FDR = FDP_val,
    power = power_val,
    rule = rule
    # coverage = mean((test_ite_true < Pval_ite$upper) & (test_ite_true > Pval_ite$lower)),
    # CI_width = mean(Pval_ite$upper - Pval_ite$lower)
  )
}

summarize_metrics_DR_Pval <- function(pval_list, sim_list, gamma,
                                      alpha_FDR = 0.05,   
                                      FDR_method = c("None","BH","Bonferroni","mLORD"),
                                      w0 = 0.14, b0 = 0.01, decay = 0.3,
                                      xi_mod = NULL
){
  purrr::map2_dfr(
    pval_list, sim_list, 
    function(pv_tbl, sim){
      pvals <- pv_tbl$p_value
      ite_true <- sim$Y1t - sim$Y0t
      
      evaluate_metrics_pval(
        Pval_ite      = pvals,
        test_ite_true = ite_true,
        gamma         = gamma,
        alpha_FDR     = alpha_FDR,
        FDR_method = FDR_method,
        w0 = w0, b0 = b0, decay = decay,
        xi_mod = xi_mod
      )
    }, .id = "rep")
}


method_mapper <- function(label, split_ID, split_ED){
  # stopifnot(label %in% c("ID_reg","ID_rf","ID_sl","ED_reg","ED_rf","ED_sl"))
  stopifnot(label %in% c(
    "ID_reg","ID_rf","ID_sl","ID_bart","ID_bma",
    "ED_reg","ED_rf","ED_sl","ED_bart","ED_bma"
  ))
  if (startsWith(label, "ID_")) {
    split_obj <- split_ID
  } else {
    split_obj <- split_ED
  }
  
  # if (endsWith(label, "reg")) {
  #   fits <- list(nuis = "reg", train="reg")
  # } else if (endsWith(label, "rf")) {
  #   fits <- list(nuis = "rf",  train="rf")
  # } else {
  #   fits <- list(nuis = "sl",  train="sl")
  # }
  if (endsWith(label, "reg")) {
    fits <- list(nuis = "reg",  train = "reg")
  } else if (endsWith(label, "rf")) {
    fits <- list(nuis = "rf",   train = "rf")
  } else if (endsWith(label, "sl")) {
    fits <- list(nuis = "sl",   train = "sl")
  } else if (endsWith(label, "bart")) {
    fits <- list(nuis = "bart", train = "bart")
  } else if (endsWith(label, "bma")) {
    fits <- list(nuis = "bma",  train = "bma")
  } else {
    stop("Unknown label suffix.")
  }
  list(split = split_obj, fit_nuis = fits$nuis, fit_train = fits$train)
}

# Build split data for one case
prepare_splits_for_case <- function(case){
  load(glue("data/simdata_{case}.RData"))         # simdata: Y, Y00, Y11, Ynull
  load(glue("data/simdata_main_{case}.RData"))    # simdata_main: S, X, A, Y0, Y1, mu0, mu1
  load(glue("data/simdata_test_{case}.RData"))    # simdata_test: testing set
  load("data/setup.RData")
  list2env(setup[case, ], envir = .GlobalEnv)
  
  set.seed(seed)
  
  rct_list <- vector("list", length(simdata_main))
  ed_list  <- vector("list", length(simdata_main))
  for (i in seq_along(simdata_main)) {
    main_i <- simdata_main[[i]]
    full_i <- simdata[[i]]
    rct_list[[i]] <- subset_by_S(main = main_i, full = full_i, group = 1)
    ed_list[[i]]  <- subset_by_S(main = main_i, full = full_i, group = 0)
  }
  
  split_data_ID <- vector("list", length(rct_list))
  split_data_ED <- vector("list", length(rct_list))
  for (i in seq_along(rct_list)) {
    split_data_ID[[i]] <- list(
      split_dat1 = split_data(data = rct_list[[i]], a = 1, use_ed = FALSE,
                              ed = NULL, borrow_trt = FALSE,
                              split_train = 0.5, split_nuisance = 0.2, method = "DR"),
      split_dat0 = split_data(data = rct_list[[i]], a = 0, use_ed = FALSE,
                              ed = NULL, borrow_trt = FALSE,
                              split_train = 0.5, split_nuisance = 0.2, method = "DR")
    )
    split_data_ED[[i]] <- list(
      split_dat1 = split_data(data = rct_list[[i]], a = 1, use_ed = TRUE,
                              ed = ed_list[[i]], borrow_trt = borrow_trt,
                              split_train = 0.5, split_nuisance = 0.2, method = "DR"),
      split_dat0 = split_data(data = rct_list[[i]], a = 0, use_ed = TRUE,
                              ed = ed_list[[i]], borrow_trt = borrow_trt,
                              split_train = 0.5, split_nuisance = 0.2, method = "DR")
    )
  }
  list(split_ID = split_data_ID, split_ED = split_data_ED, sim_test = simdata_test)
}

# Run all 6 methods for a given scenario
run_methods_for_scenario <- function(methods, split_ID, split_ED, sim_test,
                                     cs_type, cv_nuisance, cv_model, cap,
                                     alpha = alpha_overall,
                                     parallel = FALSE,
                                     workers = NULL,
                                     strategy = c("multisession","multicore")){
  pkgs_needed <- c("SuperLearner","nnls","ranger","grf","quantreg",
                   "caret","dplyr","tibble","purrr","glue", "BAS", "BART")
  if (parallel) {
    if (is.null(workers)) workers <- future::availableCores()
    old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)
    if (strategy == "multisession") {
      future::plan(multisession, workers = workers)
    } else {
      future::plan(multicore, workers = workers)   # Linux only
    }
    mapper2  <- furrr::future_map2
    map_opts <- furrr::furrr_options(seed = 2025, packages = pkgs_needed)
  } else {
    mapper2  <- purrr::map2
    map_opts <- NULL
  }
  
  Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
             OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
             BLIS_NUM_THREADS="1")
  #grf::set_num_threads(1)
  
  out_list <- vector("list", length(methods))
  names(out_list) <- methods
  
  for (m in methods){
    cfg <- method_mapper(m, split_ID, split_ED)
    tic(glue("{m} | cs={cs_type} | cv_nuis={cv_nuisance} | cv_pred={cv_model} | cap={cap} | parallel={parallel}"))
    
    res <- mapper2(
      cfg$split, sim_test,
      ~ get_Pval_DR(.x, .y,
                    alpha = alpha,
                    fit_nuisance = cfg$fit_nuis,
                    fit_train    = cfg$fit_train,
                    CV_nusiance  = cv_nuisance,
                    CV_model = cv_model,
                    conformal_type = cs_type,
                    cap = cap),
      .options = map_opts
    )
    toc()
    out_list[[m]] <- res
  }
  out_list
}

# Summarize to a table
summarize_methods <- function(results_list, sim_test, gamma = alpha_overall,
                              alpha_FDR,
                              FDR_methods = c("None","BH","Bonferroni","mLORD"),
                              w0 = 0.14, b0 = 0.01, decay = 0.3,
                              xi_mod = NULL
){
  # results_list: named list of length 6 methods, each = list of reps
  FDR_methods <- match.arg(FDR_methods, c("None","BH","Bonferroni","mLORD"), several.ok = TRUE)
  meth_names <- names(results_list)
  tidy <- purrr::map_dfr(meth_names, function(m) {
    purrr::map_dfr(FDR_methods, function(fm) {
      tbl <- summarize_metrics_DR_Pval(
        pval_list   = results_list[[m]],
        sim_list    = sim_test,
        gamma       = gamma,
        alpha_FDR   = alpha_FDR,
        FDR_method  = fm,
        w0 = w0, b0 = b0, decay = decay,
        xi_mod = xi_mod
      )
      dplyr::mutate(tbl, method = m, .before = 1)
    })
  })
  # wide summary (means per method)
  wide <- tidy %>%
    dplyr::group_by(method, rule) %>%
    dplyr::summarise(across(where(is.numeric), ~mean(.x, na.rm=TRUE)), .groups = "drop") %>%
    arrange(method,rule)
  list(tidy = tidy, wide = wide)
}

generate_xi_mod <- function(K, alpha, b0, decay){
  base <- 1 / ((1 + log(1:K)) * (1:K)^(decay))    # shape
  C <- (alpha / b0) / sum(base * (1 + log(1:K))) # scale to satisfy constraint
  xi <- C * base
  xi
}
run_modified_lord <- function(pvals, alpha = 0.05, w0 = 0.005, b0 = 0.045, decay,
                              xi = NULL) {
  n <- length(pvals)

  # sum xi_i (1 + log i) <= alpha/b0
  if (is.null(xi)) {
    xi <- generate_xi_mod(n, alpha, b0, decay)
  } else {
    if (length(xi) < n) {
      xi <- generate_xi_mod(n, alpha, b0, decay)
    }else{
      xi <- xi
    }
  }

  W <- numeric(n + 1)
  W[1] <- w0
  R <- integer(n)
  alpha_i <- numeric(n)
  tau <- 0  # last discovery time (0 = before start)

  for (i in seq_len(n)) {
    W_tau <- if (tau == 0) w0 else W[tau + 1] # wealth at last discovery time

    alpha_i[i] <- xi[i] * W_tau
    alpha_i[i] <- min(alpha_i[i], W[i])
    R[i] <- as.integer(pvals[i] <= alpha_i[i])
    W[i + 1] <- W[i] - alpha_i[i] + b0 * R[i]

    if (R[i] == 1) tau <- i
  }

  tibble(i = 1:n, pval = pvals, alpha_i = alpha_i, R_i = R, W_i = W[-1])
}
