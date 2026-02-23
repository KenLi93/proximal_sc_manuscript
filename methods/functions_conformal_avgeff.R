
###########################################################
##### kendrick.li@stjude.org                          #####
##### 02/23/2021                                      #####
###########################################################

#############################################################
## GMM method for the CI of the average effect
## Negative control approach to synthetic control, without covariates
## Input data is a list that contains the following elements:
## X: treatment indicator with length T
## Z: a T-by-N matrix of proxies
## W: a T-by-N matrix of control units
## Y: a T-by-N matrix of outcomes
#############################################################
NC_nocov_gmm <- function(data, q = 10, alpha = 0.1) {
  n.W <- ncol(data$W)
  S1 <- with(data, cbind(X, W))
  S2 <- with(data, cbind(X, Z * (1 - X)))
  Y <- data$Y
  spsz <- length(data$X)
  # obj_func <- function(tt) {
  #   c(t(Y - S1 %*% tt) %*% S2 %*% t(S2) %*% (Y - S1 %*% tt))
  # }
  # optim_rslt <- optim(par = rep(0, 1 + n.W), fn = obj_func, method = "BFGS")
  # theta_est <- optim_rslt$par
  theta_est <- c(solve(t(S1) %*% S2 %*% t(S2) %*% S1) %*% 
                   t(S1) %*% S2 %*% t(S2) %*% Y)
  bg <- c(Y - S1 %*% theta_est) * S2
  bG <- - solve(t(S2) %*% S1 / spsz)
  
  # "sandwich" variance estimate and Newey-West HAC variance
  hacOmega <- Omega <- t(bg) %*% bg / spsz
  for(i in 1:q){
    Omega_i <- t(bg[-(1:i),]) %*% bg[1:(spsz-i),] / spsz
    hacOmega <- hacOmega + (1 - i/(q+1))*(Omega_i + t(Omega_i))
  }
  Sigma <- bG %*% Omega %*% t(bG)
  hacSigma <- bG %*% hacOmega %*% t(bG)
  
  VAR <- Sigma/spsz; HACVAR <- hacSigma/spsz
  CI_lb <- theta_est[1] - qnorm(1 - alpha / 2) * sqrt(VAR[1, 1])
  CI_ub <- theta_est[1] + qnorm(1 - alpha / 2) * sqrt(VAR[1, 1])
  
  HAC_CI_lb <- theta_est[1] - qnorm(1 - alpha / 2) * sqrt(HACVAR[1, 1])
  HAC_CI_ub <- theta_est[1] + qnorm(1 - alpha / 2) * sqrt(HACVAR[1, 1])
  return(c(CI_lb = CI_lb, CI_ub = CI_ub, 
           HAC_CI_lb = HAC_CI_lb, HAC_CI_ub = HAC_CI_ub))
}


#############################################################
## GMM method for the CI of the average effect
## OLS to synthetic control, without covariates
## Input data is a list that contains the following elements:
## X: treatment indicator with length T
## Z: a T-by-N matrix of proxies
## W: a T-by-N matrix of control units
## Y: a T-by-N matrix of outcomes
#############################################################
OLS_nocov_gmm <- function(data, q = 10, alpha = 0.1) {
  n.W <- ncol(data$W)
  S1 <- with(data, cbind(X, W, Z))
  S2 <- with(data, cbind(X, W * (1 - X), Z * (1 - X)))
  Y <- data$Y
  spsz <- length(data$X)
  # obj_func <- function(tt) {
  #   c(t(Y - S1 %*% tt) %*% S2 %*% t(S2) %*% (Y - S1 %*% tt))
  # }
  # optim_rslt <- optim(par = rep(0, 1 + n.W), fn = obj_func, method = "BFGS")
  # theta_est <- optim_rslt$par
  theta_est <- c(solve(t(S1) %*% S2 %*% t(S2) %*% S1) %*% 
                   t(S1) %*% S2 %*% t(S2) %*% Y)
  bg <- c(Y - S1 %*% theta_est) * S2
  bG <- - solve(t(S2) %*% S1 / spsz)
  
  # "sandwich" variance estimate and Newey-West HAC variance
  hacOmega <- Omega <- t(bg) %*% bg / spsz
  for(i in 1:q){
    Omega_i <- t(bg[-(1:i),]) %*% bg[1:(spsz-i),] / spsz
    hacOmega <- hacOmega + (1 - i/(q+1))*(Omega_i + t(Omega_i))
  }
  Sigma <- bG %*% Omega %*% t(bG)
  hacSigma <- bG %*% hacOmega %*% t(bG)
  
  VAR <- Sigma/spsz; HACVAR <- hacSigma/spsz
  CI_lb <- theta_est[1] - qnorm(1 - alpha / 2) * sqrt(VAR[1, 1])
  CI_ub <- theta_est[1] + qnorm(1 - alpha / 2) * sqrt(VAR[1, 1])
  
  HAC_CI_lb <- theta_est[1] - qnorm(1 - alpha / 2) * sqrt(HACVAR[1, 1])
  HAC_CI_ub <- theta_est[1] + qnorm(1 - alpha / 2) * sqrt(HACVAR[1, 1])
  return(c(CI_lb = CI_lb, CI_ub = CI_ub, 
           HAC_CI_lb = HAC_CI_lb, HAC_CI_ub = HAC_CI_ub))
}

##############################################################
## SCPI method for prediction interval of the average effect
## Input data is a list that contains the following elements:
## X: treatment indicator with length T
## Z: a T-by-N matrix of proxies
## W: a T-by-N matrix of control units
## Y: a T-by-N matrix of outcomes
#############################################################
scpi_sc_int <- function(data, t0, method = "ols", alpha = 0.1) {
  Z <- data$Z
  W <- data$W
  Y <- data$Y
  tt <- length(Y)
  
  
  Y_dat <- data.frame(val = Y, time = 1:tt, id = 1)
  W_list <- apply(W, 2, function(ww) {
    data.frame(val = ww, time = 1:tt)
  })
  for (nW in 1:ncol(W)) {
    W_list[[nW]]$id <- 1 + nW
  }
  W_dat <- dplyr::bind_rows(W_list)
  
  Z_list <- apply(Z, 2, function(zz) {
    data.frame(val = zz, time = 1:tt)
  })
  for (nZ in 1:ncol(Z)) {
    Z_list[[nZ]]$id <- 1 + ncol(W) + nZ
  }
  Z_dat <- dplyr::bind_rows(Z_list)
  
  new_dat <- rbind(Y_dat, W_dat, Z_dat)
  
  sc_df <- scdata(df = new_dat, 
                  id.var = "id", 
                  time.var = "time",
                  outcome.var = "val", 
                  period.pre = 1:t0, 
                  period.post = (t0 + 1):tt, 
                  unit.tr = 1,
                  unit.co = 2:(1 + ncol(Z) + ncol(W)))
  if (method == "simplex") {
    res_pi <- scpi.ate(data = sc_df, 
                       w.constr = list(name = "simplex"),
                       u.sigma = "HC1",
                       e.method = "gaussian",
                       e.alpha = alpha / 2,
                       u.alpha = alpha / 2,
                       V.mat = diag(t0) / t0)
  } else if (method == "ols") {
    res_pi <- scpi.ate(data = sc_df, 
                       w.constr = list(name = "ols"),
                       u.sigma = "HC1",
                       e.method = "gaussian",
                       e.alpha = alpha / 2,
                       u.alpha = alpha / 2,
                       V.mat = diag(t0) / t0)
  }
  eff_pi <- mean(Y[(t0 + 1):tt]) - as.numeric(res_pi$inference.results$CI.all.gaussian[1, 2:1])
  return(eff_pi)
}

##############################################################
## SCPI method for prediction interval of the average effect
## The proximal inference approach
## Input data is a list that contains the following elements:
## X: treatment indicator with length T
## Z: a T-by-N matrix of proxies
## W: a T-by-N matrix of control units
## Y: a T-by-N matrix of outcomes
#############################################################
scpi_nc_int <- function(data, t0, method = "ols", alpha = 0.1) {
  Z <- data$Z
  W <- data$W
  Y <- data$Y
  tt <- length(Y)
  
  
  Y_dat <- data.frame(val = Y, time = 1:tt, id = 1)
  W_list <- apply(W, 2, function(ww) {
    data.frame(val = ww, time = 1:tt)
  })
  for (nW in 1:ncol(W)) {
    W_list[[nW]]$id <- 1 + nW
  }
  W_dat <- dplyr::bind_rows(W_list)
  
  new_dat <- rbind(Y_dat, W_dat)
  
  sc_df <- scdata(df = new_dat, 
                  id.var = "id", 
                  time.var = "time",
                  outcome.var = "val", 
                  period.pre = 1:t0, 
                  period.post = (t0 + 1):tt, 
                  unit.tr = 1,
                  unit.co = 2:(1 + ncol(W)))
  if (method == "simplex") {
    res_pi <- scpi.ate(data = sc_df, 
                       w.constr = list(name = "simplex"),
                       u.sigma = "HC1",
                       e.method = "gaussian",
                       e.alpha = alpha / 2,
                       u.alpha = alpha / 2,
                       V.mat = Z[1:t0, ] %*% t(Z[1:t0, ]) / t0)
  } else if (method == "ols") {
    res_pi <- scpi.ate(data = sc_df, 
                       w.constr = list(name = "ols"),
                       u.sigma = "HC1",
                       e.method = "gaussian",
                       e.alpha = alpha / 2,
                       u.alpha = alpha / 2,
                       V.mat = Z[1:t0, ] %*% t(Z[1:t0, ]) / t0)
  }
  eff_pi <- mean(Y[(t0 + 1):tt]) - as.numeric(res_pi$inference.results$CI.all.gaussian[1, 2:1])
  return(eff_pi)
}




