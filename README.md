## Code for replication of the analysis in "Theory for Identification and Inference with Synthetic Controls: A Proximal Causal Inference Framework"

This repository contains R scripts to replicate the simulation and data application in the manuscript "Theory for Identification and Inference with Synthetic Controls: A Proximal Causal Inference Framework". The simulation was performed on high-performing computing (HPC) system using LSF.

Below are the descriptions of the directories and files:

- simulation: 
  - run_sim_linear_est.R, run_sim_conformal_PI.R, run_sim_avgeff_PI.R: R scripts for generating the data and performing the simulation studies;
  - functions_linear.R, functions_conformal.R, functions_conformal_avgeff.R, scpi_ate.R: methods to perform effect estimation or conformal inference of the effect size estimates;
  - plot_LM.R, plot_TV.R, tab_lm_N30.R: summarize and visualize the simulation results.
- application:
  - german_data_functions.R: utility functions for preprocessing the German reunification data;
  - functions_conformal.R, functions_conformal_avgeff.R, scpi_ate.R: methods to perform effect estimation or conformal inference of the effect size estimates;
  - german_data_app.R: the code for data analysis using the German reunification data.
 
For any questions, please contact Kendrick Li (kendrick.li@stjude.org).
