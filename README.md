## Code for replication of the analysis in "Theory for Identification and Inference with Synthetic Controls: A Proximal Causal Inference Framework"

This repository contains R scripts to replicate the simulation and data application in the manuscript "Theory for Identification and Inference with Synthetic Controls: A Proximal Causal Inference Framework". The simulation was performed on high-performing computing (HPC) system using LSF.

Below are the descriptions of the directories and files:

# simulation: 
  - run_sim_linear_est.R, run_sim_conformal_PI.R, run_sim_avgeff_PI.R: R scripts for generating the data and performing the simulation studies. 
  
  - plot_LM.R, plot_TV.R, tab_lm_N30.R, summarize_ate_pi.R, summarize_pointwise_pi.R: summarize and visualize the simulation results.
    - Figure 2 and Supplemental Figures S.1 - S.5 are generated using results from run_sim_linear_est.R, followed by code in plot_LM.R;
    - Figures S.6 and S.7 are generated using results from run_sim_linear_est.R, followed by code in plot_TV.R;
    - Table S.1 are generated using results from run_sim_linear_est.R, followed by code in tab_lm_N30.R;
    - Tables 1 and S.2-S.3 is generated using results from  run_sim_conformal_PI.R, followed by summarize_pointwise_pi.R;
    - Table 2 is generated using results from run_sim_avgeff_PI.R, followed by code in summarize_ate_pi.R.
# application:
  - german_data_functions.R: utility functions for preprocessing the German reunification data;
  - functions_conformal.R, functions_conformal_avgeff.R, scpi_ate.R: methods to perform effect estimation or conformal inference of the effect size estimates;
  - german_data_app.R: the code for data analysis using the German reunification data.
    - Figures 3 and S.9 are generated using code in german_data_app.R.
# methods:
  - functions_linear.R, functions_conformal.R, functions_conformal_avgeff.R, scpi_ate.R: methods to perform effect estimation or conformal inference of the effect size estimates.


For any questions, please contact Kendrick Li (kendrick.li@stjude.org).
