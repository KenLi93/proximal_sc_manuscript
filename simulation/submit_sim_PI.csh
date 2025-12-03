#!/bin/tcsh

#BSUB -P 1117
#BSUB -J sim_SC_PI[11-120]
#BSUB -q gpu
#BSUB -n 1
#BSUB -gpu "num=1"
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=64GB]"
#BSUB -u kli@stjude.org
#BSUB -B
#BSUB -N
#BSUB -e ./Rout/sim_PI_err_%J_%I.err

module load R/4.3.2
Rscript run_sim_conformal_PI.R 
