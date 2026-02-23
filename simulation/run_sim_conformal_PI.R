###########################################################
##### Pointwise predictive inference under the linear #####
##### interactive fixed effects model using           #####
##### synthetic control methods                       #####
##### (1) Permutation inference method                #####
##### (2) SCPI method                                 #####
##### kendrick.li@stjude.org                          #####
##### 02/23/2021                                      #####
###########################################################
rm(list=ls())

library(Synth)
library(sandwich)
library(Rsolnp)
library(dplyr)
library(parallel)
library(scpi)
source("functions_conformal.R")
param_grid <- expand.grid(n.rep = 500,
                          #t0 = c(20, 80, 140, 200),
                          t0 = c(30),
                          #n.units = c(1 + 2, 1 + 10, 1 + 20),
                          n.units = c(5, 7, 11),
                          dist.epsilon = c("iid"),
                          dist.lambda = c("stationary", "nonstationary"),
                          #dist.lambda = "nonstationary",
                          U.setting = c("constrained", "unconstrained"),
                          addcov = c(FALSE),
                          alpha = c(0.1),
                          batch = 1:10)



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
run.one <- function(seed, n.units, t0, mysd, theta, true.beta, addcov = F,
                    dist.epsilon, dist.lambda, U.setting, epsilonARMAcor = 0.1,
                    alpha = 0.1, grid_lb = -25, grid_ub = 25, output = "both"){
  set.seed(seed)
  .env <- environment()
  
  n.ctrl <- n.units - 1 ## one treated unit
  n.Z <- floor(n.ctrl / 2) ## use half of the control units as proxies
  n.W <- n.ctrl - n.Z ## number of donors
  ind.trt <- 1 ## index of treated units
  ind.ctrl <- setdiff(1:n.units, ind.trt)
  ind.Z <- 1:n.Z ## index of proxies among the control units
  
  ## indicator of treatment
  X <<- c(rep(0, t0), rep(1, t0))
  
  ### gen lambda and U
  U <- generate.U(n.units, 2 * t0, U.setting)
  
  lambda.t <- gen.lambda(n.units, 2 * t0, dist.lambda)
  
  ### generate covariates
  if(addcov == T){
    C <- replicate(n = n.units, expr = rnorm(2 * t0, mean = 0, sd = 1))
  }
  
  ### generate outcome
  Y.allunits = NULL
  
  for(i in 1:n.units){
    ### first generate stationary weakly dependent error term epsilon
    if(dist.epsilon=="iid"){
      epsilon = rnorm(2 * t0, mean = 0, sd = mysd)
    }else if(dist.epsilon=="AR"){
      epsilon = as.numeric(arima.sim(n = 2 * t0, list(ar = epsilonARMAcor), innov = rnorm(2 * t0)))
    }else{
      print('unrecognized dist.epsilon')
      break
    }
    ### then generate outcome for all units
    if(addcov == T){
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
  
  
  data <- list(Y = Y[1:(t0 + 1)], Z = Z[1:(t0 + 1), ], W = W[1:(t0 + 1), ])
  data.full <- list(Y = Y, Z = Z, W = W)
  crude_CI <- ConformalPointwiseCI(data = data, t0 = t0, grid = seq(grid_lb, grid_ub, 0.1),
                                   addcov = addcov, output = "both", 
                                   method = "NC", alpha = alpha)
  ## fine search
  NC_CI_lb <- 
    ConformalPointwiseCI(data = data, t0 = t0, 
                         grid = seq(crude_CI[1] - 0.09, crude_CI[1], by = 0.01),
                         addcov = addcov, output = "lb", 
                         method = "NC", alpha = alpha)
  
  NC_CI_ub <- 
    ConformalPointwiseCI(data = data, t0 = t0, 
                         grid = seq(crude_CI[2], crude_CI[2] + 0.09, by = 0.01),
                         addcov = addcov, output = "ub", 
                         method = "NC", alpha = alpha)
  
  ## NC constrained
  crude_CI <- ConformalPointwiseCI(data = data, t0 = t0, grid = seq(grid_lb, grid_ub, 0.1),
                                   addcov = addcov, output = "both", 
                                   method = "NC-constrained", alpha = alpha)
  ## fine search
  NC_constrained_CI_lb <- 
    ConformalPointwiseCI(data = data, t0 = t0, 
                         grid = seq(crude_CI[1] - 0.09, crude_CI[1], by = 0.01),
                         addcov = addcov, output = "lb", 
                         method = "NC-constrained", alpha = alpha)
  
  NC_constrained_CI_ub <- 
    ConformalPointwiseCI(data = data, t0 = t0, 
                         grid = seq(crude_CI[2], crude_CI[2] + 0.09, by = 0.01),
                         addcov = addcov, output = "ub", 
                         method = "NC-constrained", alpha = alpha)
  ### OLS method
  ## crude search
  crude_CI <- ConformalPointwiseCI(data = data, t0 = t0, grid = seq(grid_lb, grid_ub, 0.1),
                                   addcov = addcov, output = "both", 
                                   method = "OLS", alpha = alpha)
  ## fine search
  OLS_CI_lb <- 
    ConformalPointwiseCI(data = data, t0 = t0, 
                         grid = seq(crude_CI[1] - 0.09, crude_CI[1], by = 0.01),
                         addcov = addcov, output = "lb", 
                         method = "OLS", alpha = alpha)
  OLS_CI_ub <- 
    ConformalPointwiseCI(data = data, t0 = t0, 
                         grid = seq(crude_CI[2], crude_CI[2] + 0.09, by = 0.01),
                         addcov = addcov, output = "ub", 
                         method = "OLS", alpha = alpha)
  ### SC method
  ## crude search
  crude_CI <- ConformalPointwiseCI(data = data, t0 = t0, grid = seq(grid_lb, grid_ub, 0.1),
                                   addcov = addcov, output = "both", 
                                   method = "SC", alpha = alpha)
  ## fine search
  SC_CI_lb <- 
    ConformalPointwiseCI(data = data, t0 = t0, 
                         grid = seq(crude_CI[1] - 0.09, crude_CI[1], by = 0.01),
                         addcov = addcov, output = "lb", 
                         method = "SC", alpha = alpha)
  SC_CI_ub <- 
    ConformalPointwiseCI(data = data, t0 = t0, 
                         grid = seq(crude_CI[2], crude_CI[2] + 0.09, by = 0.01),
                         addcov = addcov, output = "ub", 
                         method = "SC", alpha = alpha)
  ## SCPI 
  scpi_sc_constrained <- scpi_sc_int(data = data.full, t0 = t0, method = "simplex", alpha = alpha)
  scpi_sc_unconstrained <- scpi_sc_int(data = data.full, t0 = t0, method = "ols", alpha = alpha)
  
  scpi_nc_constrained <- scpi_nc_int(data = data.full, t0 = t0, method = "simplex", alpha = alpha)
  scpi_nc_unconstrained <- scpi_nc_int(data = data.full, t0 = t0, method = "ols", alpha = alpha)
  
  
  return(c(NC_lb = NC_CI_lb, NC_ub = NC_CI_ub,
           NC_constrained_lb = NC_constrained_CI_lb, NC_constrained_ub = NC_constrained_CI_ub,
           OLS_lb = OLS_CI_lb, OLS_ub = OLS_CI_ub,
           SC_lb = SC_CI_lb, SC_ub = SC_CI_ub,
           scpi_sc_constrained_lb = scpi_sc_constrained[1], scpi_sc_constrained_ub = scpi_sc_constrained[2],
           scpi_sc_unconstrained_lb = scpi_sc_unconstrained[1], scpi_sc_unconstrained_ub = scpi_sc_unconstrained[2],
           scpi_nc_constrained_lb = scpi_nc_constrained[1], scpi_nc_constrained_ub = scpi_nc_constrained[2],
           scpi_nc_unconstrained_lb = scpi_nc_unconstrained[1], scpi_nc_unconstrained_ub = scpi_nc_unconstrained[2]))
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
alpha <- param_grid[job_id, "alpha"]
batch <- param_grid[job_id, "batch"]




mysd <- 1.5

theta <- 1
true.beta <- 2

epsilonARMAcor <- 0.1


## save results
myseeds <- (batch - 1) * n.rep + 1:n.rep

rslt.all <- mclapply(myseeds,
                   function(seed) {
                     print(seed)
                     
                     run.one(seed = seed, n.units = n.units, t0 = t0,
                             mysd = mysd, theta = theta, true.beta = true.beta,
                             addcov = addcov, dist.epsilon = dist.epsilon,
                             dist.lambda = dist.lambda, U.setting = U.setting,
                             epsilonARMAcor = epsilonARMAcor,
                             alpha = alpha)
                     
                   }, mc.cores = detectCores())

saveRDS(rslt.all, file = paste0("results_PI/LM",
                                "_epsilon_", dist.epsilon,
                                "_U_", U.setting,
                                "_lambda_", dist.lambda,
                                "_cov_", addcov,
                                "_mysd_", mysd,
                                "_t0_", t0,
                                "_nunits_", n.units,
                                "_nrep_", n.rep,
                                "_alpha_", alpha,
                                "_batch_", batch,
                                ".rds"))
