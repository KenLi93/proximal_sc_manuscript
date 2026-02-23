###########################################################
##### Summary table with T0 = 30                      #####
##### Setting: Linear interactive fixed effects model #####
##### Synthetic control methods                       #####
##### (1) Constrained OLS (Abadie)                    #####
##### (2) Unconstrained OLS                           #####
##### (3) Proximal inference                          #####
##### kendrick.li@stjude.org                          #####
##### 02/23/2026                                      #####
###########################################################
rm(list=ls())

rm(list=ls())

library(dplyr)
library(xtable)
########################################################################
#################################### bias plot #########################
########################################################################
# fixed.lambda <- F; 
true.beta <- 2; n.rep <- 500
t0.all <- c(30)
#n.units.all <- c(1 + 4, 1 + 10, 1 + 20);
n.units.all <- c(1 + 4, 1 + 6, 1 + 10);
xbound <- 0.005
# n.units.all=c(1+10,1+20,1+30);


# U_is_fixed <- T; 
myrange <- 200

# if (mysd == 0.5) { myrange <- 100 }

param_grid <- rbind(expand.grid(dist.epsilon = c("iid"),
                                dist.lambda = c("stationary", "nonstationary"),
                                U.setting = c("constrained", "unconstrained"),
                                addcov = c(TRUE, FALSE)),
                    expand.grid(dist.epsilon = c("AR"),
                                dist.lambda = c("stationary", "nonstationary"),
                                U.setting = c("constrained", "unconstrained"),
                                addcov = c(FALSE))) %>% unique()




## Reading all datasets and organizing into a single data frame
rslt.all.combine <- NULL

for (ii in 1:nrow(param_grid)) {      
  addcov <- param_grid[ii, "addcov"]
  dist.epsilon <- param_grid[ii, "dist.epsilon"]
  dist.lambda <- param_grid[ii, "dist.lambda"]
  U.setting <- param_grid[ii, "U.setting"]
  
  if (dist.lambda == "nonstationary") {
    mysd <- 1.5
  } else {
    mysd <- 1.5
  }
    # for(arg3 in c(1,2,4,6,7)[1]){
  
  for (n.units in n.units.all) {
    for (t0 in t0.all) {
      for(batch in 1:10){
        tmp_list <- try(readRDS(file = paste0("results_linear_est/LM_v5",
                                              "_epsilon_", dist.epsilon,
                                              "_U_", U.setting,
                                              "_lambda_", dist.lambda,
                                              "_cov_", addcov,
                                              "_mysd_", mysd,
                                              "_t0_", t0,
                                              "_nunits_", n.units,
                                              "_nrep_", n.rep,
                                              "_batch_", batch,
                                              ".rds")),
                        silent=T)
        
        if (class(tmp_list)[[1]] != "try-error"){
          tmp_list <- tmp_list[unlist(sapply(tmp_list, 
                                             function(x) class(x) != "try-error"))]
          tmp <- data.frame(t0 = t0,
                            addcov = addcov,
                            dist.lambda = dist.lambda,
                            U.setting = U.setting,
                            n.units = n.units,
                            SC_est = sapply(tmp_list, function(x) x["SC_est"]),
                            SC_OLS_est = sapply(tmp_list, function(x) x["SC_OLS_est"]),
                            NC_SC_est = sapply(tmp_list, function(x) x["NC_SC_est"]),
                            NC_SC_constrained_est = sapply(tmp_list, function(x) x["NC_SC_constrained_est"]),
                            SC_OLS_se = sapply(tmp_list, function(x) x["SC_OLS_se"]),
                            SC_OLS_se_hac = sapply(tmp_list, function(x) x["SC_OLS_se_hac"]),
                            NC_SC_se = sapply(tmp_list, function(x) x["NC_SC_se"]),
                            NC_SC_se_hac = sapply(tmp_list, function(x) x["NC_SC_se_hac"]),
                            SC_est2 = sapply(tmp_list, function(x) x["SC_est2"]),
                            SC_OLS_est2 = sapply(tmp_list, function(x) x["SC_OLS_est2"]),
                            NC_SC_est2 = sapply(tmp_list, function(x) x["NC_SC_est2"]),
                            NC_SC_constrained_est2 = sapply(tmp_list, function(x) x["NC_SC_constrained_est2"]),
                            SC_OLS_se2 = sapply(tmp_list, function(x) x["SC_OLS_se2"]),
                            SC_OLS_se_hac2 = sapply(tmp_list, function(x) x["SC_OLS_se_hac2"]),
                            NC_SC_se2 = sapply(tmp_list, function(x) x["NC_SC_se2"]),
                            NC_SC_se_hac2 = sapply(tmp_list, function(x) x["NC_SC_se_hac2"]))
          
          rslt.all.combine <- rbind(rslt.all.combine, tmp)
        }
      } 
    }
  }
}

## For each scenario, obtain the bias, median absolute error, and SD of
## different estimators

rslt.all.summary.nocov <- rslt.all.combine %>%
  filter(addcov == F) %>%
  group_by(dist.lambda, U.setting, n.units) %>%
  summarise(OLS_bias = mean(SC_OLS_est - 2),
            OLS_mae = median(SC_OLS_est - 2),
            OLS_sd = sd(SC_OLS_est),
            SC_bias = mean(SC_est - 2),
            SC_mae = median(SC_est - 2),
            SC_sd = sd(SC_est),
            PI_bias = mean(NC_SC_est - 2),
            PI_mae = median(NC_SC_est - 2),
            PI_sd = sd(NC_SC_est),
            cPI_bias = mean(NC_SC_constrained_est - 2),
            cPI_mae = median(NC_SC_constrained_est - 2),
            cPI_sd = sd(NC_SC_constrained_est)) %>%
  mutate(across(OLS_bias:cPI_sd, ~ format(round(., 2), nsmall = 2))) %>%
  arrange(dist.lambda, U.setting, n.units)

rslt.nocov.xtab <- xtable(rslt.all.summary.nocov, digits = rep(0, ncol(rslt.all.summary.nocov) + 1))

print(rslt.nocov.xtab, file = "rslt_nocov_xtable.txt", include.rownames = F)



# ########################################################################
# #################################### 95%CI ####################################
# ########################################################################
# rm(list=ls())
# 
# qq <- 1 - (1 - 0.95) / 2 ##95%CI
# 
# 
# param_grid <- expand.grid(n.rep = 200,
#                           t0 = c(50, 100, 200),
#                           n.units = c(1 + 2, 1 + 10, 1 + 20),
#                           addcov = c(T, F),
#                           dist.epsilon = c("iid", "AR"))
# true.beta <- 2
# rslt.all.combine <- NULL
# for (ii in 1:nrow(param_grid)) {
#   n.rep <- param_grid[ii, "n.rep"]
#   t0 <- param_grid[ii, "t0"]
#   n.units <- param_grid[ii, "n.units"]
#   dist.epsilon <- param_grid[ii, "dist.epsilon"]
#   addcov <- param_grid[ii, "addcov"]
#   
#   for (batch in 1:10) {
#     if (dist.epsilon == "iid") {
#       rslt <- readRDS(file = paste0("results/LM_",
#                                     "epsilon_", dist.epsilon,
#                                     "cov_", addcov,
#                                     "t0_", t0,
#                                     "nunits_", n.units,
#                                     "nrep", n.rep,
#                                     "batch", batch,
#                                     ".rds")) %>% 
#         as.data.frame %>%
#         transmute(SC.cover = as.numeric(SC.ate - qnorm(qq) * SC.se < true.beta &
#                                           SC.ate + qnorm(qq) * SC.se > true.beta ),
#                   NC.cover = as.numeric(NC.ate2 - qnorm(qq) * NC.se2 < true.beta &
#                                           NC.ate2 + qnorm(qq) * NC.se2 > true.beta ),
#                   CORR = dist.epsilon,
#                   addcov = addcov,
#                   nunits = n.units,
#                   t0 = t0)
#       rslt.all.combine <- rbind(rslt.all.combine, rslt)
#     } else {
#       rslt <- readRDS(file = paste0("results/LM_",
#                                     "epsilon_", dist.epsilon,
#                                     "cov_", addcov,
#                                     "t0_", t0,
#                                     "nunits_", n.units,
#                                     "nrep", n.rep,
#                                     "batch", batch,
#                                     ".rds")) %>% 
#         as.data.frame %>%
#         transmute(SC.cover = as.numeric(SC.ate - qnorm(qq) * SC.se < true.beta &
#                                           SC.ate + qnorm(qq) * SC.se > true.beta ),
#                   NC.cover = as.numeric(NC.ate2 - qnorm(qq) * NC.se2.HAC < true.beta &
#                                           NC.ate2 + qnorm(qq) * NC.se2.HAC > true.beta ),
#                   CORR = dist.epsilon,
#                   addcov = addcov,
#                   nunits = n.units,
#                   t0 = t0)
#       rslt.all.combine <- rbind(rslt.all.combine, rslt)
#     } 
#   }
# }
# 
# rslt.all.summary <- rslt.all.combine %>% 
#   group_by(CORR, nunits, addcov, t0) %>%
#   summarize(SC.coverage = mean(SC.cover), NC.coverage = mean(NC.cover))
# 
# write.csv(rslt.all.summary, file = "tab_LM_coverage.csv")
# 
# 
 