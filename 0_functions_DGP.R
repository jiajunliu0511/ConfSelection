# to make ITE model be more flexible
ite_linear <- function(X, sd_noise_spec, beta_ite=c(2,3), c0=0.6){#0.6(Dec11);0.7;c0=0.8
  X1 <- X[,1]
  X2 <- X[,2]
  as.vector(X %*% beta_ite) + c0 + X1 * X2
}
ite_nonlinear <- function(X, sd_noise_spec) {
  X1 <- X[,1]
  X2 <- X[,2]
  
  # base <- 3.5 * (X1 * X2 + exp(X2 - 1))
  base <- 3.5 * X1 * X2 + 4*exp(1.5*(X2 - 1))
  shift <- dplyr::case_when(
    sd_noise_spec == "const:0.3"   ~ -1.1,
    sd_noise_spec == "const:0.6"   ~ -0.8,
    sd_noise_spec == "hetero:ITE"  ~ -1.0,
    TRUE                           ~ -1.0   # default / safeguard
  )
  as.vector(base + shift)
}
get_ite_fun <- function(ite_type) {
  switch(ite_type,
         "linear"    = ite_linear,
         "nonlinear" = ite_nonlinear,
         stop("Unknown ite_type: ", ite_type))
}

# noise structure setup
noise_sd_from_spec <- function(ITE, sd_noise_spec = "const:0.5") {
  if (startsWith(sd_noise_spec, "const:")) {
    k_str <- substr(sd_noise_spec, 7, nchar(sd_noise_spec))
    k <- as.numeric(k_str)
    if (!is.na(k)) {
      return(rep(k, length(ITE)))
    } else {
      stop("Invalid constant value in sd_noise_spec: ", sd_noise_spec)
    }
  }
  if (sd_noise_spec == "hetero:ITE") {
    sd_vec <- 0.3 * (0.25*ITE^2*(abs(ITE)<2) + 0.5*abs(ITE)*(abs(ITE)>=1))
    bump <- 0.5 * exp(-(ITE / 0.9)^2) # 0.1*exp... is too easy
    return(pmax(sd_vec+bump, 0.2)) #
    # sd_vec <- 0.3 * (0.25*ITE^2*(abs(ITE)<2) + 0.5*abs(ITE)*(abs(ITE)>=1))
    # return(pmax(sd_vec, 1e-6))
  }
  stop("Unknown sd_noise_spec: ", sd_noise_spec)
}

# generate the data
generate_main <- function(n,
                          etaj = 0.3, # covariate shift
                          borrow_trt = TRUE,
                          main_only = TRUE, # FALSE: generate ED
                          ite_fun,
                          sd_noise_spec = NULL,
                          seed = NULL) {    
  if (!is.null(seed)) set.seed(seed)
  
  if (main_only) {
    # ---- Generate testing data & RCT----
    X <- cbind(runif(n, -1, 1), runif(n, -1, 1))
    colnames(X) <- c("X1","X2")
    X1 <- X[,1]
    X2 <- X[,2]
    
    mu0 <- 0.5 + X1 + X2
    ITE <- ite_fun(X, sd_noise_spec = sd_noise_spec)
    mu1  <- mu0 + ITE
    
    sd_noise <- noise_sd_from_spec(ITE, sd_noise_spec)
    Y0 <- mu0 + rnorm(n, 0, sd_noise)
    Y1 <- mu1 + rnorm(n, 0, sd_noise)
    R2 <- mean(c(var(mu0) / var(Y0), var(mu1) / var(Y1)))
    A <- rbinom(n, 1, 0.5)
    S <- rep(1, n)
    
    return(list(
      X = X, S = S, A = A,
      mu0 = mu0, mu1 = mu1, ITE = ITE,
      Y0 = Y0, Y1 = Y1,
      R2 = R2,
      sd_noise_spec = sd_noise_spec, 
      sd_noise = sd_noise
    ))
  }else{
    # ---- Generate EC ----
    N_super <- max(3 * n, 10000)
    get_one_superpop <- function() {
      X <- cbind(runif(N_super, -1, 1), runif(N_super, -1, 1))
      colnames(X) <- c("X1","X2")
      X1 <- X[, 1]; X2 <- X[, 2]
      
      eta  <- rep(etaj, 2)
      Xeta <- as.vector(X %*% eta)
      target_pS1 <- (N_super - n) / N_super
      
      froot <- function(a0) mean(plogis(a0 + Xeta)) - target_pS1
      a0    <- uniroot(froot, c(-100, 100))$root
      
      piS <- plogis(a0 + Xeta)
      S   <- rbinom(N_super, 1, piS)
      
      list(X = X, X1 = X1, X2 = X2, S = S)
    }
    success <- FALSE
    while (!success) {
      super <- get_one_superpop()
      idx_S0 <- which(super$S == 0)
      if (length(idx_S0) >= n) {
        success <- TRUE
      }
    }
    keep <- sample(idx_S0, n)
    
    X  <- super$X[keep, , drop = FALSE]
    X1 <- X[, 1]; X2 <- X[, 2]
    S  <- rep(0, n)
    A <- if (borrow_trt) rbinom(n, 1, 0.4) else rep(0, n)
    mu0 <- 0.5 + X1 + X2
    ITE <- ite_fun(X, sd_noise_spec = sd_noise_spec)
    mu1 <- mu0 + ITE
    
    sd_noise <- noise_sd_from_spec(ITE, sd_noise_spec)
    Y0 <- mu0 + rnorm(n, 0, sd_noise)
    Y1 <- mu1 + rnorm(n, 0, sd_noise)
    R2 <- mean(c(var(mu0) / var(Y0), var(mu1) / var(Y1)))
    
    return(list(X = X, S = S, A = A, 
                mu0 = mu0, mu1 = mu1, ITE = ITE, 
                Y0 = Y0, Y1 = Y1, 
                R2 = R2, 
                sd_noise_spec = sd_noise_spec, sd_noise = sd_noise))
  }
}

# consider when there is outcome drift
generate_bias_from_main <- function(main_item,
                                borrow_trt = TRUE,
                                mag_bias = 0.2,
                                prop_unbias = 0.5,
                                sd_noise_spec = NULL,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  mu0 <- main_item$mu0
  mu1 <- main_item$mu1
  S   <- main_item$S
  A   <- main_item$A
  
  N <- length(S)
  
  s0_idx <- which(S == 0)
  if (mag_bias == 0) {
    mu00 <- mu0
    mu11 <- if (borrow_trt) mu1 else rep(NA_real_, N)
    id_unbias <- s0_idx
  } else {
    unbias_n <- round(length(s0_idx) * prop_unbias)
    id_unbias <- if (unbias_n > 0) sample(s0_idx, unbias_n) else integer(0)
    
    mu00 <- mu0 - mag_bias
    # mu11 <- if (borrow_trt) mu1 - mag_bias else rep(NA_real_, N)
    mu11 <- if (borrow_trt) mu1 else rep(NA_real_, N)
    
    if (length(id_unbias) > 0) {
      mu00[id_unbias] <- mu0[id_unbias]
      if (borrow_trt) mu11[id_unbias] <- mu1[id_unbias]
    }
  }
  ITE_hb <- mu11 - mu00
  sd_noise <- noise_sd_from_spec(ITE_hb, sd_noise_spec)
  Y00 <- mu00 + rnorm(N, 0, sd_noise)
  Y11 <- if (borrow_trt) mu11 + rnorm(N, 0, sd_noise) else rep(NA_real_, N)
  
  
  # Observed Y and the all-control Ynull
  Y <- S*A*main_item$Y1 + S*(1 - A)*main_item$Y0 + (1 - S)*(1 - A)*Y00 + (1 - S)*A*Y11
  Ynull <- S*main_item$Y0 + (1 - S)*Y00
  
  list(mu00 = mu00, mu11 = mu11, ITE_hb = ITE_hb, Y00 = Y00, Y11 = Y11, Y = Y, Ynull = Ynull, id_unbias = id_unbias, sd_noise_spec = sd_noise_spec, sd_noise = sd_noise)
}

