#!/bin/bash
#SBATCH --job-name=Congo
#SBATCH --partition=fuchs
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=8:00:00

# Activate conda
source /home/fuchs/fias/knguyen/.bashrc
conda activate kinh
cd /scratch/fuchs/fias/knguyen/MultistageSurv/

srun R CMD BATCH --no-save --no-restore ./code/run.r

