###########################################################
##### Setting: Linear interactive fixed effects model #####
##### Synthetic control methods                       #####
##### (1) Constrained OLS (Abadie)                    #####
##### (2) Unconstrained OLS                           #####
##### (3) Proximal inference                          #####
##### kendrick.li@stjude.org                          #####
##### 02/23/2026                                      #####
###########################################################
rm(list=ls())


library(sandwich)
library(Rsolnp)
library(dplyr)
library(parallel)
library(scpi)
source("functions_linear.R")
param_grid <- rbind(expand.grid(n.rep = 500,
                                t0 = c(30, 80, 140, 200),
                                #n.units = c(1 + 2, 1 + 10, 1 + 20),
                                n.units = c(5, 7, 11),
                                dist.epsilon = c("iid"),
                                dist.lambda = c("stationary", "nonstationary"),
                                #dist.lambda = "nonstationary",
                                U.setting = c("constrained", "unconstrained"),
                                addcov = c(TRUE, FALSE),
                                batch = 1:10),
                    expand.grid(n.rep = 500,
                                t0 = c(30, 80, 100, 200),
                                # n.units = c(1 + 2, 1 + 10, 1 + 20),
                                n.units = c(5, 7, 11),
                                dist.epsilon = c("AR"),
                                dist.lambda = c("stationary", "nonstationary"),
                                # dist.lambda = "nonstationary",
                                U.setting = c("constrained", "unconstrained"),
                                addcov = c(FALSE),
                                batch = 1:10)) %>% unique()


## latent factors 
gen.lambda <- function(n.units, t, dist.lambda){
  n.lambda <- floor((n.units - 1) / 2) ## number of latent factors
  if (dist.lambda == "stationary") {
    mylambda <- cbind(replicate(n = n.lambda, rnorm(t, 0.5, sd = 0.5)))
  } else {
    mylambda <- cbind(replicate(n = n.lambda, 0.5 * log(1:t) + rnorm(t, sd = 0.5)))
  }
  return(mylambda)
}

## factor loading
generate.U <- function(n.units, t, U.setting){
  n.Z <- floor((n.units-1) / 2) ## one treated unit, half of the control units are proxies
  n.W <- n.units - 1 - n.Z
  if (U.setting == "unconstrained") {
    U_0 <- 1.5
    U.mat <- matrix(0, ncol = n.W, nrow = n.W)
    
    ele_diag <- c(2, rep(1, n.W - 1))
    # ele_diag <- 1:n.W
    diag(U.mat) <- sum(ele_diag) / ele_diag
  } else if (U.setting == "constrained") {
    U_0 <- 1
    U.mat <- matrix(0, ncol = n.W, nrow = n.W)
    
    ele_diag <- c(2, rep(1, n.W - 1))
    # ele_diag <- 1:n.W
    diag(U.mat) <- sum(ele_diag) / ele_diag
  }
  U <- rbind(U_0, U.mat, U.mat)
  return(U)
}

## implement one run of simulation
run.one <- function(seed, n.units, t, t0, mysd, theta, true.beta, addcov = F,
                    dist.epsilon, dist.lambda, U.setting, epsilonARMAcor = 0.1){
  set.seed(seed)
  .env <- environment()
  
  n.ctrl <- n.units - 1 ## one treated unit
  n.Z <- floor(n.ctrl / 2) ## use half of the control units as proxies
  n.W <- n.ctrl - n.Z ## number of donors
  ind.trt <- 1 ## index of treated units
  ind.ctrl <- setdiff(1:n.units, ind.trt)
  ind.Z <- 1:n.Z ## index of proxies among the control units
  
  ## indicator of treatment
  X <<- c(rep(0, t0), rep(1, t - t0))
  
  ### gen lambda and U
  U <- generate.U(n.units, t, U.setting)
  
  lambda.t <- gen.lambda(n.units, t, dist.lambda)
  
  ### generate covariates
  if(addcov == T){
    C <- replicate(n = n.units, expr = rnorm(t, mean = 0, sd = 1))
  }
  
  ### generate outcome
  Y.allunits = NULL
  
  for(i in 1:n.units){
    ### first generate stationary weakly dependent error term epsilon
    if(dist.epsilon=="iid"){
      epsilon = rnorm(t, mean = 0, sd = mysd)
    }else if(dist.epsilon=="AR"){
      epsilon = as.numeric(arima.sim(n = t, list(ar = epsilonARMAcor), innov = rnorm(t)))
    }else{
      print('unrecognized dist.epsilon')
      break
    }
    ### then generate outcome for all units
    if(addcov==T){
      Yi = c(C[,i] * theta + U[i, ] %*% t(lambda.t) + epsilon)
    }else{
      Yi = c(U[i, ] %*% t(lambda.t) + epsilon)
    }
    Y.allunits = cbind(Y.allunits, Yi)
  }
  colnames(Y.allunits)=1:(n.ctrl+1)
  
  ##
  Y <- Y.allunits[,ind.trt] + true.beta * X
  V <- Y.allunits[,ind.ctrl]
  Z <- cbind(Y.allunits[,ind.ctrl[ind.Z]])
  W <- cbind(Y.allunits[,ind.ctrl[-ind.Z]])
  U.Y <- cbind(U[ind.trt,])
  U.Z <- t(U[ind.ctrl[ind.Z],])
  U.W <- t(U[ind.ctrl[-ind.Z],])
  if(addcov == T){
    C.Y <- cbind(C[,ind.trt])
    C.Z <- cbind(C[,ind.ctrl[ind.Z]])
    C.W <- cbind(C[,ind.ctrl[-ind.Z]])
  }
  #### define vcov assumption on epsilon
  if(dist.epsilon == "iid"){
    vcov.epsilon <- "iid"
  } else {
    vcov.epsilon <- "HAC"
  }
  
  
  #### SC methods
  if(addcov == T){
    data <- list(X = X, Y = Y, W = W, Z = Z, C.Y = C.Y, C.W = C.W, C.Z = C.Z)
    
    SC_rslt <- SC_cov(data)
    SC_OLS_rslt <- SC_OLS_cov(data)
    NC_SC_rslt <- NC_SC_cov(data)
    NC_SC_constrained_rslt <- NC_SC_constrained_cov(data)
    
    SC_rslt2 <- SC_nocov(data)
    SC_OLS_rslt2 <- SC_OLS_nocov(data)
    NC_SC_rslt2 <- NC_SC_nocov(data)
    NC_SC_constrained_rslt2 <- NC_SC_constrained_nocov(data)
  } else {
    data <- list(X = X, Y = Y, W = W, Z = Z)
    
    SC_rslt <- SC_nocov(data)
    SC_OLS_rslt <- SC_OLS_nocov(data)
    NC_SC_rslt <- NC_SC_nocov(data)
    NC_SC_constrained_rslt <- NC_SC_constrained_nocov(data)
  }
  
  
  
  SC_est <- SC_rslt$est
  SC_OLS_est <- SC_OLS_rslt$est
  NC_SC_est <- NC_SC_rslt$est
  NC_SC_constrained_est <- NC_SC_constrained_rslt$est
  
  SC_OLS_se <- SC_OLS_rslt$se
  SC_OLS_se_hac <- SC_OLS_rslt$se_hac
  NC_SC_se <- NC_SC_rslt$se
  NC_SC_se_hac <- NC_SC_rslt$se_hac
  if (addcov == TRUE) {
    SC_est2 <- SC_rslt2$est
    SC_OLS_est2 <- SC_OLS_rslt2$est
    NC_SC_est2 <- NC_SC_rslt2$est
    NC_SC_constrained_est2 <- NC_SC_constrained_rslt2$est
    
    SC_OLS_se2 <- SC_OLS_rslt$se
    SC_OLS_se_hac2 <- SC_OLS_rslt$se_hac
    NC_SC_se2 <- NC_SC_rslt$se
    NC_SC_se_hac2 <- NC_SC_rslt$se_hac
  } else {
    SC_est2 <- SC_OLS_est2 <- NC_SC_est2 <- NC_SC_constrained_est2 <- 
      SC_OLS_se2 <- SC_OLS_se_hac2 <- NC_SC_se2 <- NC_SC_se_hac2 <- NA
  }
  return(c(SC_est = SC_est, 
           SC_OLS_est = SC_OLS_est,
           NC_SC_est = NC_SC_est, 
           NC_SC_constrained_est = NC_SC_constrained_est,
           SC_OLS_se = SC_OLS_se,
           SC_OLS_se_hac = SC_OLS_se_hac,
           NC_SC_se = NC_SC_se,
           NC_SC_se_hac = NC_SC_se_hac,
           SC_est2 = SC_est2, 
           SC_OLS_est2 = SC_OLS_est2,
           NC_SC_est2 = NC_SC_est2, 
           NC_SC_constrained_est2 = NC_SC_constrained_est2,
           SC_OLS_se2 = SC_OLS_se2,
           SC_OLS_se_hac2 = SC_OLS_se_hac2,
           NC_SC_se2 = NC_SC_se2,
           NC_SC_se_hac2 = NC_SC_se_hac2))
}


#job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
job_id <- as.numeric(Sys.getenv('LSB_JOBINDEX'))
print(job_id)
n.rep <- param_grid[job_id, "n.rep"]
t0 <- param_grid[job_id, "t0"]
n.units <- param_grid[job_id, "n.units"]
dist.epsilon <- param_grid[job_id, "dist.epsilon"]
dist.lambda <- param_grid[job_id, "dist.lambda"]
U.setting <- param_grid[job_id, "U.setting"]
addcov <- param_grid[job_id, "addcov"]
batch <- param_grid[job_id, "batch"]




mysd <- 1.5

theta <- 1
true.beta <- 2

epsilonARMAcor <- 0.1

t_star <- t0
t <- t_star + t0

## save results
myseeds <- (batch - 1) * n.rep + 1:n.rep

rslt.all <- mclapply(myseeds,
                     function(seed) {
                       print(seed)
                       
                       run.one(seed = seed, n.units = n.units, t = t, t0 = t0,
                               mysd = mysd, theta = theta, true.beta = true.beta,
                               addcov = addcov, dist.epsilon = dist.epsilon,
                               dist.lambda = dist.lambda, U.setting = U.setting,
                               epsilonARMAcor = epsilonARMAcor)
                       
                     }, mc.cores = detectCores())

saveRDS(rslt.all, file = paste0("results_linear_est/LM_v5",
                                "_epsilon_", dist.epsilon,
                                "_U_", U.setting,
                                "_lambda_", dist.lambda,
                                "_cov_", addcov,
                                "_mysd_", mysd,
                                "_t0_", t0,
                                "_nunits_", n.units,
                                "_nrep_", n.rep,
                                "_batch_", batch,
                                ".rds"))
