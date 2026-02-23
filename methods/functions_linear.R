#############################################################
## classical synthetic control method, with covariates
## Input data is a list that contains the following elements:
## X: treatment indicator with length T
## Z: a T-by-N matrix of proxies
## W: a T-by-N matrix of control units
## Y: a T-by-N matrix of outcomes
## C.Y: a T-by-N matrix of covariates for the outcome units
## C.W: a T-by-N matrix of covariates for the control units
## C.Z: a T-by-N matrix of covariates for the proxy units
#############################################################

SC_cov <- function(data) {
  t0 <- sum(1 - data$X)
  tt <- length(data$X)
  
  Z <- data$Z
  W <- data$W
  Y <- data$Y
  C.Y <- data$C.Y
  C.W <- data$C.W
  C.Z <- data$C.Z
  X <- data$X
  tt <- length(Y)
  
  
  Y_dat <- data.frame(val = Y, cov = C.Y, time = 1:tt, id = 1)
  W_list <- vector("list", ncol(W))
  for (nW in 1:ncol(W)) {
    W_list[[nW]] <- data.frame(val = W[, nW], cov = C.W[, nW], time = 1:tt, id = 1 + nW)
  }
  W_dat <- dplyr::bind_rows(W_list)
  
  Z_list <- vector("list", ncol(Z))
  for (nZ in 1:ncol(Z)) {
    Z_list[[nZ]] <- data.frame(val = Z[, nZ], cov = C.Z[, nZ], time = 1:tt, id = 1 + ncol(W) + nZ)
  }
  Z_dat <- dplyr::bind_rows(Z_list)
  
  new_dat <- rbind(Y_dat, W_dat, Z_dat)
  
  sc_df <- scdata(df = new_dat, 
                  id.var = "id", 
                  time.var = "time",
                  outcome.var = "val", 
                  period.pre = 1:t0, 
                  period.post = (t0 + 1):tt, 
                  cov.adj = list("cov" = c("constant", "trend")),
                  unit.tr = 1,
                  unit.co = 2:(1 + ncol(Z) + ncol(W)))
  
  res_sc <- scest(data = sc_df, 
                  w.constr = list(name = "simplex"))
  
  ww_unordered <- res_sc$est.results$w
  ww_ordered <- ww_unordered[order(as.numeric(stringr::str_extract(names(ww_unordered), "\\d+\\z")))]
  
  SC.ate.constrained <-
    mean(Y[(t0+1):tt] - c(cbind(W, Z)[(t0+1):tt, ] %*% ww_ordered))
  
  
  return(list(est = SC.ate.constrained,
              SC_wts = ww_ordered))
}

#############################################################
## classical synthetic control method, without covariates
## Input data is a list that contains the following elements:
## X: treatment indicator with length T
## Z: a T-by-N matrix of proxies
## W: a T-by-N matrix of control units
## Y: a T-by-N matrix of outcomes
#############################################################
SC_nocov <- function(data) {
  Z <- data$Z
  W <- data$W
  Y <- data$Y
  X <- data$X
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
                  period.pre = which(X == 0), 
                  period.post = which(X == 1), 
                  unit.tr = 1,
                  unit.co = 2:(1 + ncol(Z) + ncol(W)))
  
  res_sc <- scest(data = sc_df, 
                  w.constr = list(name = "simplex"))
  
  ww_unordered <- res_sc$est.results$w
  ww_ordered <- ww_unordered[order(as.numeric(stringr::str_extract(names(ww_unordered), "\\d+\\z")))]
  
  ate <- mean(Y[which(X == 1)] - c(cbind(W, Z)[which(X == 1),] %*% ww_ordered))
  return(list(est = ate,
              SC_wts = c(ww_ordered)))
}


#############################################################
## OLS method for synthetic control with covariates
## Input data is a list that contains the following elements:
## X: treatment indicator with length T
## Z: a T-by-N matrix of proxies
## W: a T-by-N matrix of control units
## Y: a T-by-N matrix of outcomes
## C.Y: a T-by-N matrix of covariates for the outcome units
## C.W: a T-by-N matrix of covariates for the control units
## C.Z: a T-by-N matrix of covariates for the proxy units
#############################################################
SC_OLS_cov <- function(data, q = 10) {
  S1 <- with(data, cbind(X, W, Z, C.Y, C.Z, C.W))
  S2 <- with(data, cbind(X, W * (1 - X), Z * (1 - X),
                         C.Y * (1 - X), C.Z * (1 - X), C.W * (1 - X)))
  Y <- data$Y
  spsz <- length(data$X)
  
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
  
  return(list(est = theta_est[1], se = sqrt(VAR[1, 1]), 
              se_hac = sqrt(HACVAR[1, 1])))
}

#############################################################
## OLS method for SC, without covariates
## Input data is a list that contains the following elements:
## X: treatment indicator with length T
## Z: a T-by-N matrix of proxies
## W: a T-by-N matrix of control units
## Y: a T-by-N matrix of outcomes
#############################################################
SC_OLS_nocov <- function(data, q = 10) {
  S1 <- with(data, cbind(X, W, Z))
  S2 <- with(data, cbind(X, W * (1 - X), Z * (1 - X)))
  Y <- data$Y
  spsz <- length(data$X)

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
  
  return(list(est = theta_est[1], se = sqrt(VAR[1, 1]), 
              se_hac = sqrt(HACVAR[1, 1])))
}

#############################################################
## Negative control approach to synthetic control, with covariates
## Input data is a list that contains the following elements:
## X: treatment indicator with length T
## Z: a T-by-N matrix of proxies
## W: a T-by-N matrix of control units
## Y: a T-by-N matrix of outcomes
## C.Y: a T-by-N matrix of covariates for the outcome units
## C.W: a T-by-N matrix of covariates for the control units
## C.Z: a T-by-N matrix of covariates for the proxy units
#############################################################

NC_SC_cov <- function(data, q = 10) {
  n.W <- ncol(data$W)
  S1 <- with(data, cbind(X, W, C.Y, C.W))
  S2 <- with(data, cbind(X, Z * (1 - X), C.Y * (1 - X), C.W * (1 - X)))
  Y <- data$Y
  spsz <- length(data$X)
  # obj_func <- function(tt) {
  #   c(t(Y - S1 %*% tt) %*% S2 %*% t(S2) %*% (Y - S1 %*% tt))
  # }
  # optim_rslt <- optim(par = rep(0, 2 + 2 * n.W), fn = obj_func, method = "BFGS")
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
  
  return(list(est = theta_est[1], 
              SC_wts = theta_est[2:(1 + n.W)],
              coef_C = theta_est[-c(1:(1 + n.W))],
              se = sqrt(VAR[1, 1]), 
              se_hac = sqrt(HACVAR[1, 1])))
}

#############################################################
## Negative control approach to synthetic control, with covariates
## SC weights are constrained in a simplex
## Input data is a list that contains the following elements:
## X: treatment indicator with length T
## Z: a T-by-N matrix of proxies
## W: a T-by-N matrix of control units
## Y: a T-by-N matrix of outcomes
## C.Y: a T-by-N matrix of covariates for the outcome units
## C.W: a T-by-N matrix of covariates for the control units
## C.Z: a T-by-N matrix of covariates for the proxy units
#############################################################

NC_SC_constrained_cov <- function(data) {
  n.W <- ncol(data$W)
  
  S1 <- with(data, cbind(X, W, C.Y, C.W))
  S2 <- with(data, cbind(X, Z * (1 - X), C.Y * (1 - X), C.W * (1 - X)))
  Y <- data$Y
  spsz <- length(data$X)
  
  ## objective function for optimization
  if (n.W > 1) { ## have more than one control units
    eval_func <- function(theta) {
      as.numeric(t(Y - S1 %*% theta) %*% S2 %*% t(S2) %*% (Y - S1 %*% theta))
    }
    
    ## SC weights sum up to one
    eval_constraint_eq <- function(theta) {
      sum(theta[2:(1 + n.W)]) - 1
    } 
    
    eval_constraint_ineq <- function(theta) {
      theta[2:(1 + n.W)]
    }
    
    # local_opts <- list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15)
    # opts <- list( "algorithm"= "NLOPT_GN_ISRES",
    #               "xtol_rel"= 1.0e-15,
    #               "maxeval"= 160000,
    #               "local_opts" = local_opts,
    #               "print_level" = 0 )
    # 
    # using unconstrained results as initial value
    res_unconstrained <- NC_SC_cov(data)
    SC_wts_init <- pmax(res_unconstrained$SC_wts, 0.001)
    theta_init <- c(res_unconstrained$est,
                    SC_wts_init / sum(SC_wts_init),
                    res_unconstrained$coef_C)
    theta_init[2:(1 + n.W)] <- theta_init[2:(1 + n.W)] / 
      sum(theta_init[2:(1 + n.W)])
    
    
    opt_result <- solnp(pars = theta_init,
                        fun = eval_func,
                        eqfun = eval_constraint_eq,
                        eqB = 0,
                        ineqfun = eval_constraint_ineq,
                        ineqLB = rep(0, n.W),
                        ineqUB = rep(1, n.W),
                        control = list(trace = FALSE))
    
    # opt_result <- nloptr(x = theta_init,
    #                      eval_f = eval_func,
    #                      # eval_grad_f = eval_grad_func,
    #                      lb = theta_lb,
    #                      ub = theta_ub,
    #                      eval_g_eq = eval_constraint_eq,
    #                      # eval_jac_g_eq = eval_constants_eq_grad,
    #                      opts = opts)
    
    return(list(est = opt_result$pars[1], 
                SC_wts = opt_result$pars[2:(1 + n.W)],
                coef_C = opt_result$pars[-(1:(1 + n.W))]))
    
    
    
  } else {
    eval_func <- function(par) {
      est <- par[1]
      cov_eff <- par[-1]
      theta <- c(est, 1, cov_eff)
      as.numeric(t(Y - S1 %*% theta) %*% S2 %*% t(S2) %*% (Y - S1 %*% theta))
    }
    
    res_unconstrained <- NC_SC_cov(data)
    
    est_init <- res_unconstrained$est
    reg_coef_init <- res_unconstrained$coef_C
    opt_result <- optim(par = c(est_init, reg_coef_init),
                        fn = eval_func)
    return(list(est = opt_result$par[1], 
                SC_wts = 1,
                coef_C = opt_result$par[-1]))
    
  }
  
}

#############################################################
## Negative control approach to synthetic control, without covariates
## Input data is a list that contains the following elements:
## X: treatment indicator with length T
## Z: a T-by-N matrix of proxies
## W: a T-by-N matrix of control units
## Y: a T-by-N matrix of outcomes
#############################################################
NC_SC_nocov <- function(data, q = 10) {
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
  
  return(list(est = theta_est[1],
              SC_wts = theta_est[-1],
              se = sqrt(VAR[1, 1]), 
              se_hac = sqrt(HACVAR[1, 1])))
}


#############################################################
## Negative control approach to synthetic control, without covariates
## SC weights are constrained in a simplex
## Input data is a list that contains the following elements:
## X: treatment indicator with length T
## Z: a T-by-N matrix of proxies
## W: a T-by-N matrix of control units
## Y: a T-by-N matrix of outcomes
#############################################################

NC_SC_constrained_nocov <- function(data, use_posttrt_data = T) {
  n.W <- ncol(data$W)
  
  S1 <- with(data, cbind(X, W))
  S2 <- with(data, cbind(X, Z * (1 - X)))
  Y <- data$Y
  
  ## objective function for optimization
  if (n.W > 1) {
    eval_func <- function(theta) {
      as.numeric(t(Y - S1 %*% theta) %*% S2 %*% t(S2) %*% (Y - S1 %*% theta))
    } 
    
    ## gradient function
    eval_grad_func <- function(theta) {
      -2 * t(S1) %*% S2 %*% t(S2) %*% Y + 
        2 * t(S1) %*% S2 %*% t(S2) %*% S1 %*% theta
    }
    
    ## SC weights sum up to one
    eval_constraint_eq <- function(theta) {
      sum(theta[2:(1 + n.W)]) - 1
    } 
    
    eval_constraint_ineq <- function(theta) {
      theta[2:(1 + n.W)]
    }
    # eval_constants_eq_grad <- function(theta) {
    #   c(0, rep(1, n.W))
    # }
    
    # local_opts <- list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15)
    # opts <- list( "algorithm"= "NLOPT_GN_ISRES",
    #               "xtol_rel"= 1.0e-15,
    #               "maxeval"= 160000,
    #               "local_opts" = local_opts,
    #               "print_level" = 0 )
    res_unconstrained <- NC_SC_nocov(data)
    
    SC_wts_init <- pmax(res_unconstrained$SC_wts, 0.001)
    theta_init <- c(res_unconstrained$est, 
                    SC_wts_init / sum(SC_wts_init))
    
    # opt_result <- nloptr(x = theta_init,
    #                      eval_f = eval_func,
    #                      # eval_grad_f = eval_grad_func,
    #                      lb = theta_lb,
    #                      ub = theta_ub,
    #                      eval_g_eq = eval_constraint_eq,
    #                      # eval_jac_g_eq = eval_constants_eq_grad,
    #                      opts = opts)
    opt_result <- solnp(pars = theta_init,
                        fun = eval_func,
                        eqfun = eval_constraint_eq,
                        eqB = 0,
                        ineqfun = eval_constraint_ineq,
                        ineqLB = rep(0, n.W),
                        ineqUB = rep(1, n.W),
                        control = list(trace = FALSE))
    
    return(list(est = opt_result$pars[1], SC_wts = opt_result$pars[-1]))
  } else {
    eval_func <- function(par) {
      theta <- c(par, 1)
      as.numeric(t(Y - S1 %*% theta) %*% S2 %*% t(S2) %*% (Y - S1 %*% theta))
    }
    
    opt_result <- optimize(f = eval_func, interval = c(-1000, 1000))
    
    return(list(est = opt_result$minimum, SC_wts = 1))
  }
}
