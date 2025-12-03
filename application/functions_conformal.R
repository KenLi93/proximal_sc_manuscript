
## estimating equation to estimate the SC weights
NC_nocov <- function(data, t0){
  Y <- cbind(data$Y)[1:t0]
  Z <- cbind(data$Z)[1:t0,]
  W <- cbind(data$W)[1:t0,]
  
  NC_est <- try(c(solve(t(W) %*% Z %*% t(Z) %*% W) %*% t(W) %*% Z %*% t(Z) %*% Y), silent = T)
  
  if (class(NC_est[1]) == "error") {
    NC_est_init <- c(MASS::ginv(t(W) %*% Z %*% t(Z) %*% W) %*% t(W) %*% Z %*% t(Z) %*% Y)
    est_func <- function(ww) {
      uu <- Z * c(Y - W %*% ww)
      return(mean(colMeans(uu) ^ 2))
    }
    NC_est <- optim(par = NC_est_init, fn = est_func)$par
  }
  return(NC_est)
}


## estimating equation to estimate the SC weights
OLS_nocov <- function(data, t0){
  Y <- cbind(data$Y)[1:t0];
  V <- cbind(data$W, data$Z)[1:t0, ]
  
  OLS_est <- try(c(solve(t(V) %*% V) %*% t(V) %*% Y), silent = T)
  if (class(OLS_est[1]) == "error") {
    OLS_est_init <- c(MASS::ginv(t(V) %*% V) %*% t(V) %*% Y)
    est_func <- function(ww) {
      uu <- V * c(Y - V %*% ww)
      return(mean(colMeans(uu) ^ 2))
    }
    OLS_est <- optim(par = OLS_est_init, fn = est_func)$par
  }
  return(OLS_est)
}

SC_nocov <- function(data, t0) {
  Y <- data$Y
  tt <- length(Y)
  Z <- matrix(data$Z, nrow = tt)
  W <- matrix(data$W, nrow = tt)
  
  
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
  
  res_sc <- scest(data = sc_df, 
                  w.constr = list(name = "simplex"),
                  V.mat = diag(t0) / t0)
  
  ww_unordered <- res_sc$est.results$w
  ww_ordered <- ww_unordered[order(as.numeric(stringr::str_extract(names(ww_unordered), "\\d+\\z")))]
  
  return(ww_ordered)
}




NC_constrained_nocov <- function(data, t0) {
  
  ## objective function for optimization
  Y <- data$Y
  tt <- length(Y)
  Z <- matrix(data$Z, nrow = tt)
  W <- matrix(data$W, nrow = tt)
  
  
  
  
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
  
  res_sc <- scest(data = sc_df, 
                  w.constr = list(name = "simplex"),
                  V.mat = Z[1:t0, ] %*% t(Z[1:t0, ]) / t0)
  
  ww_unordered <- res_sc$est.results$w
  ww_ordered <- ww_unordered[order(as.numeric(stringr::str_extract(names(ww_unordered), "\\d+\\z")))]
  
  return(ww_ordered)
  
}

## Conformal pointwise confidence intervals


ConformalPointwiseCI <- 
  function(data, t0, t1, grid_list, output = "both", alpha = 0.1, method = "NC-constrained") {
    X <- cbind(data$X); Y <- cbind(data$Y)
    Z <- cbind(data$Z); W <<- cbind(data$W)
    V <- cbind(data$W, data$Z)
    X1 <- c(rep(0, t0), 1)
    
    lb <- ub <- 1:t1 * NA
    
    for (tt in 1:t1){
      grid <- grid_list[[tt]]
      pval <- NA * grid
      for (bb in seq_along(grid)) {
        Y1 <- cbind(c(Y[1:t0], Y[t0 + tt] - grid[bb]))
        Z1 <- rbind(Z[1:t0, ], Z[t0 + tt,])
        W1 <- rbind(W[1:t0, ], W[t0 + tt,])
        V1 <- cbind(W1, Z1)
        dat.h0 <- list(Y = Y1, Z = Z1, W = W1)
        if (method == "NC") {
          ww <- try(NC_nocov(dat.h0), silent = T)
          resid <- c(Y1  - W1 %*% ww)
        } else if (method == "OLS") {
          ww <- try(OLS_nocov(dat.h0), silent = T)
          resid <- c(Y1  - V1 %*% ww)
        } else if (method == "SC") {
          ww <- try(SC_nocov(dat.h0, t0), silent = T)
          resid <- c(Y1  - V1 %*% ww)
        } else if (method == "NC-constrained") {
          ww <- try(NC_constrained_nocov(dat.h0, t0), silent = T)
          resid <- c(Y1  - W1 %*% ww)
        }
        if (!inherits(ww[1], "try-error")) {
          
          resid_t <- abs(resid[t0 + 1])
          pval[bb] <- 1 - mean(resid_t > abs(resid))
        } else {
          pval[bb] <- NA
        }
        
        
      }
      
      # find upper bound
      if (output == "both") {
        #  SC_lb <- grid[min(which(SC.pval > alpha))]
        #  SC_ub <- grid[max(which(SC.pval > alpha))]
        lb_id <- max(which(pval[1:which.max(pval)] < alpha), na.rm = T) + 1
        ub_id <- which.max(pval) + min(which(pval[which.max(pval):length(pval)] < alpha), na.rm = T) - 2
        
        lb[tt] <- grid[lb_id]
        ub[tt] <- grid[ub_id]
      } else if (output == "lb") {
        lb_id <- min(which(pval >= alpha), na.rm = T) 
        lb[tt] <- grid[lb_id]
        
      } else if (output == "ub") {
        ub_id <- max(which(pval >= alpha), na.rm = T)
        ub[tt] <- grid[ub_id]
        
      }
    }
    return(list(lb = lb, ub = ub))
  }


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
    res_pi <- scpi(data = sc_df, 
                   w.constr = list(name = "simplex"),
                   u.sigma = "HC1",
                   e.method = "gaussian",
                   e.alpha = alpha / 2,
                   u.alpha = alpha / 2,
                   V.mat = diag(t0) / t0)
  } else if (method == "ols") {
    res_pi <- scpi(data = sc_df, 
                   w.constr = list(name = "ols"),
                   u.sigma = "HC1",
                   e.method = "gaussian",
                   e.alpha = alpha / 2,
                   u.alpha = alpha / 2,
                   V.mat = diag(t0) / t0)
  }
  lb <- Y[(t0 + 1):tt] - as.numeric(res_pi$inference.results$CI.all.gaussian[, 2])
  ub <- Y[(t0 + 1):tt] - as.numeric(res_pi$inference.results$CI.all.gaussian[, 1])
  return(list(lb, ub))
}

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
    res_pi <- scpi(data = sc_df, 
                   w.constr = list(name = "simplex"),
                   u.sigma = "HC1",
                   e.method = "gaussian",
                   e.alpha = alpha / 2,
                   u.alpha = alpha / 2,
                   V.mat = Z[1:t0, ] %*% t(Z[1:t0, ]) / t0)
  } else if (method == "ols") {
    res_pi <- scpi(data = sc_df, 
                   w.constr = list(name = "ols"),
                   u.sigma = "HC1",
                   e.method = "gaussian",
                   e.alpha = alpha / 2,
                   u.alpha = alpha / 2,
                   V.mat = Z[1:t0, ] %*% t(Z[1:t0, ]) / t0)
  }
  eff_pi <- Y[t0 + 1] - as.numeric(res_pi$inference.results$CI.all.gaussian[1, 2:1])
  return(eff_pi)
}




