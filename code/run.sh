#!/bin/bash
#SBATCH --job-name=msm8
#SBATCH --partition=fuchs
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --output=logs/%x_%A_%a.out
# One array task per line of the task file, each line "<CC> <sex>".
# Submitted by code/hpc_submit.sh, which sets --array and the dependency on compile.sh.
source /home/fuchs/fias/knguyen/.bashrc
conda activate kinh
cd "${SLURM_SUBMIT_DIR}"
TASKS="${1:-code/tasks_f.txt}"
read -r CC SEX < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKS")
echo "task ${SLURM_ARRAY_TASK_ID}: ${CC} ${SEX} on $(hostname)"
srun Rscript code/run.r "$CC" "$SEX"
