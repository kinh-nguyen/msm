#!/bin/bash
#SBATCH --job-name=msm8_compile
#SBATCH --partition=fuchs
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=0:30:00
#SBATCH --output=logs/%x_%j.out
# Compile msm.cpp once on a compute node. Submitted by code/hpc_submit.sh.
source /home/fuchs/fias/knguyen/.bashrc
conda activate kinh
cd "${SLURM_SUBMIT_DIR}"
Rscript code/compile.r
