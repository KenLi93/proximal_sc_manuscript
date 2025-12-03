rm(list = ls())

library(dplyr)
library(xtable)
param_grid <- expand.grid(n.rep = 500,
                          t0 = c(80, 140, 200, 1000),
                          #n.units = c(1 + 2, 1 + 10, 1 + 20),
                          n.units = c(5, 7, 11),
                          dist.epsilon = c("iid"),
                          dist.lambda = c("stationary", "nonstationary"),
                          #dist.lambda = "nonstationary",
                          U.setting = c("constrained", "unconstrained"),
                          addcov = c(FALSE),
                          alpha = c(0.1))

results_list <- vector("list", nrow(param_grid))

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
    res <- readRDS(file = paste0("results_PI/LM",
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

results_summary <- bind_rows(results_list) %>%
  group_by(t0, n.units, dist.lambda, U.setting) %>%
  summarise(NC_cover = mean(NC_lb < 2 & NC_ub > 2), 
            NC_len = mean(NC_ub - NC_lb),
            NC_constrained_cover = mean(NC_constrained_lb < 2 & NC_constrained_ub > 2),
            NC_constrained_len = mean(NC_constrained_ub - NC_constrained_lb),
            OLS_cover = mean(OLS_lb < 2 & OLS_ub > 2),
            OLS_len = mean(OLS_ub - OLS_lb),
            SC_cover = mean(SC_lb < 2 & SC_ub > 2),
            SC_len = mean(SC_ub - SC_lb),
            scpi_sc_cover = mean(scpi_sc_constrained_lb < 2 & scpi_sc_constrained_ub > 2),
            scpi_sc_len = mean(scpi_sc_constrained_ub - scpi_sc_constrained_lb),
            scpi_ols_cover = mean(scpi_sc_unconstrained_lb < 2 & scpi_sc_unconstrained_ub > 2),
            scpi_ols_len = mean(scpi_sc_unconstrained_ub - scpi_sc_unconstrained_lb),
            scpi_nc_constrained_cover = mean(scpi_nc_constrained_lb < 2 & scpi_nc_constrained_ub > 2),
            scpi_nc_constrained_len = mean(scpi_nc_constrained_ub - scpi_nc_constrained_lb),
            scpi_nc_cover = mean(scpi_nc_unconstrained_lb < 2 & scpi_nc_unconstrained_ub > 2),
            scpi_nc_len = mean(scpi_nc_unconstrained_ub - scpi_nc_unconstrained_lb)) %>%
  ungroup() %>%
  mutate(across(5:20, function(x) round(x, 3)))

main_tab <- results_summary %>%
  filter(dist.lambda == "stationary", U.setting == "constrained", t0 <= 200) %>%
  transmute(n.units, t0, 
            perm_OLS = paste0(OLS_cover * 100, "% (", round(OLS_len, 1), ")"),
            perm_SC = paste0(SC_cover * 100, "% (", round(SC_len, 1), ")"),
            perm_PI= paste0(NC_cover * 100, "% (", round(NC_len, 1), ")"),
            perm_cPI= paste0(NC_constrained_cover * 100, "% (", round(NC_constrained_len, 1), ")"),
            scpi_OLS = paste0(scpi_ols_cover * 100, "% (", round(scpi_ols_len, 1), ")"),
            scpi_SC = paste0(scpi_sc_cover * 100, "% (", round(scpi_sc_len, 1), ")"),
            scpi_PI = paste0(scpi_nc_cover * 100, "% (", round(scpi_nc_len, 1), ")"),
            scpi_cPI = paste0(scpi_nc_constrained_cover * 100, "% (", round(scpi_nc_constrained_len, 1), ")")) %>%
  arrange(n.units, t0)

main_xtab <- xtable(main_tab)
print.xtable(main_xtab, include.rownames = F)

sup1_tab <- results_summary %>%
  filter(dist.lambda == "nonstationary", U.setting == "constrained", t0 <= 200) %>%
  transmute(n.units, t0, 
            perm_OLS = paste0(OLS_cover * 100, "% (", round(OLS_len, 1), ")"),
            perm_SC = paste0(SC_cover * 100, "% (", round(SC_len, 1), ")"),
            perm_PI= paste0(NC_cover * 100, "% (", round(NC_len, 1), ")"),
            perm_cPI= paste0(NC_constrained_cover * 100, "% (", round(NC_constrained_len, 1), ")"),
            scpi_OLS = paste0(scpi_ols_cover * 100, "% (", round(scpi_ols_len, 1), ")"),
            scpi_SC = paste0(scpi_sc_cover * 100, "% (", round(scpi_sc_len, 1), ")"),
            scpi_PI = paste0(scpi_nc_cover * 100, "% (", round(scpi_nc_len, 1), ")"),
            scpi_cPI = paste0(scpi_nc_constrained_cover * 100, "% (", round(scpi_nc_constrained_len, 1), ")")) %>%
  arrange(n.units, t0)
sup1_xtab <- xtable(sup1_tab, digits = c(0, 0, 0, rep(0, 8)))
print.xtable(sup1_xtab, include.rownames = F)


sup2_tab <- results_summary %>%
  filter(dist.lambda == "stationary", U.setting == "unconstrained", t0 <= 200) %>%
  transmute(n.units, t0, 
            perm_OLS = paste0(OLS_cover * 100, "% (", round(OLS_len, 1), ")"),
            perm_SC = paste0(SC_cover * 100, "% (", round(SC_len, 1), ")"),
            perm_PI= paste0(NC_cover * 100, "% (", round(NC_len, 1), ")"),
            perm_cPI= paste0(NC_constrained_cover * 100, "% (", round(NC_constrained_len, 1), ")"),
            scpi_OLS = paste0(scpi_ols_cover * 100, "% (", round(scpi_ols_len, 1), ")"),
            scpi_SC = paste0(scpi_sc_cover * 100, "% (", round(scpi_sc_len, 1), ")"),
            scpi_PI = paste0(scpi_nc_cover * 100, "% (", round(scpi_nc_len, 1), ")"),
            scpi_cPI = paste0(scpi_nc_constrained_cover * 100, "% (", round(scpi_nc_constrained_len, 1), ")")) %>%
  arrange(n.units, t0)
sup2_xtab <- xtable(sup2_tab, digits = c(0, 0, 0, rep(0, 8)))
print.xtable(sup2_xtab, include.rownames = F)


write.csv(results_summary, file = "pointwise_pi_results.csv", row.names = F)
