#!/bin/tcsh

#BSUB -P 1117
#BSUB -J sim_linear[1-360]
#BSUB -q gpu
#BSUB -n 1
#BSUB -gpu "num=1"
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=32GB]"
#BSUB -u kli@stjude.org
#BSUB -B
#BSUB -N
#BSUB -e ./Rout/sim_PI_err_%J_%I.err
#BSUB -o ./Rout/sim_PI_err_%J_%I.out

module load R/4.3.2
Rscript run_sim_linear_est.R 
