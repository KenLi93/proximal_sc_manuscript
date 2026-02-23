###########################################################
##### Summary table                                   #####
##### Predictive inference for average effect under   #####
##### the linear interactive fixed effects model      #####
##### using synthetic control methods                 #####
##### (1) GMM method                                  #####
##### (2) SCPI method                                 #####
##### kendrick.li@stjude.org                          #####
##### 02/23/2021                                      #####
###########################################################

rm(list = ls())

library(dplyr)

param_grid <- expand.grid(n.rep = 500,
                          t0 = c(80, 140, 200),
                          #n.units = c(1 + 2, 1 + 10, 1 + 20),
                          n.units = c(5, 7, 11),
                          dist.epsilon = c("iid"),
                          dist.lambda = c("stationary", "nonstationary"),
                          #dist.lambda = "nonstationary",
                          U.setting = c("constrained", "unconstrained"),
                          addcov = c(FALSE),
                          alpha = c(0.1))

results_list <- vector("list", nrow(param_grid))

## Reading simulation results and organizing into a single data frame
for (i in 1:length(results_list)) {
  temp_df <- NULL
  n.rep <- param_grid[i, "n.rep"]
  t0 <- param_grid[i, "t0"]
  n.units <- param_grid[i, "n.units"]
  dist.epsilon <- param_grid[i, "dist.epsilon"]
  dist.lambda <- param_grid[i, "dist.lambda"]
  U.setting <- param_grid[i, "U.setting"]
  addcov <- param_grid[i, "addcov"]
  alpha <- param_grid[i, "alpha"]
  
  for (bb in 1:10) {
    res <- readRDS(file =  paste0("results_PI_ate/LM",
                                  "_epsilon_", dist.epsilon,
                                  "_U_", U.setting,
                                  "_lambda_", dist.lambda,
                                  "_cov_", addcov,
                                  "_mysd_", 1.5,
                                  "_t0_", t0,
                                  "_nunits_", n.units,
                                  "_nrep_", n.rep,
                                  "_alpha_", alpha,
                                  "_batch_", bb,
                                  ".rds"))
    
    res_df <- bind_rows(lapply(res, function(xx) {
      if (class(xx)[1] != "try-error") {
        return(xx)
      }
    }))
    
    temp_df <- rbind(temp_df, res_df)
  }
  
  results_list[[i]] <- cbind(temp_df, t0, n.units, dist.lambda, U.setting)
}

## Make the summary table for coverage probabilities and average interval lengths 
results_summary <- bind_rows(results_list) %>%
  filter(NC_gmm_lb > -20, NC_gmm_ub < 20, 
         NC_scpi_unconstrained_lb > -20, NC_scpi_unconstrained_ub < 20) %>%
  group_by(t0, n.units, dist.lambda, U.setting) %>%
  summarise(NC_gmm_cover = mean(NC_gmm_lb < 2 & NC_gmm_ub > 2, na.rm = T), 
            NC_gmm_len = mean(NC_gmm_ub - NC_gmm_lb, na.rm = T),
            OLS_gmm_cover = mean(OLS_gmm_lb < 2 & OLS_gmm_ub > 2, na.rm = T), 
            OLS_gmm_len = mean(OLS_gmm_ub - OLS_gmm_lb, na.rm = T),
            NC_scpi_constrained_cover = mean(NC_scpi_constrained_lb < 2 & NC_scpi_constrained_ub > 2, na.rm = T),
            NC_scpi_constrained_len = mean(NC_scpi_constrained_ub - NC_scpi_constrained_lb, na.rm = T),
            NC_scpi_unconstrained_cover = mean(NC_scpi_unconstrained_lb < 2 & NC_scpi_unconstrained_ub > 2, na.rm = T),
            NC_scpi_unconstrained_len = mean(NC_scpi_unconstrained_ub - NC_scpi_unconstrained_lb, na.rm = T),
            SC_scpi_cover = mean(SC_scpi_constrained_lb < 2 & SC_scpi_constrained_ub > 2, na.rm = T),
            SC_scpi_len = mean(SC_scpi_constrained_ub - SC_scpi_constrained_lb, na.rm = T),
            OLS_scpi_cover = mean(SC_scpi_unconstrained_lb < 2 & SC_scpi_unconstrained_ub > 2, na.rm = T),
            OLS_scpi_len = mean(SC_scpi_unconstrained_ub - SC_scpi_unconstrained_lb, na.rm = T)) %>%
  ungroup() %>%
  mutate(across(5:16, function(x) round(x, 3)))


main_tab <- results_summary %>%
  filter(dist.lambda == "stationary", U.setting == "constrained", t0 <= 200) %>%
  transmute(n.units, t0, 
            OLS_gmm = paste0(OLS_gmm_cover * 100, "% (", round(OLS_gmm_len, 1), ")"),
            NC_gmm = paste0(NC_gmm_cover * 100, "% (", round(NC_gmm_len, 1), ")"),
            scpi_OLS = paste0(OLS_scpi_cover * 100, "% (", round(OLS_scpi_len, 1), ")"),
            scpi_SC = paste0(SC_scpi_cover * 100, "% (", round(SC_scpi_len, 1), ")"),
            scpi_PI = paste0(NC_scpi_unconstrained_cover * 100, "% (", round(NC_scpi_unconstrained_len, 1), ")"),
            scpi_cPI = paste0(NC_scpi_constrained_cover * 100, "% (", round(NC_scpi_constrained_len, 1), ")")) %>%
  arrange(n.units, t0) %>% t

main_xtab <- xtable(main_tab)
print.xtable(main_xtab, include.rownames = T)


write.csv(results_summary, file = "ate_pi_results.csv", row.names = F)
