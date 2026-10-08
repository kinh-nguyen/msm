#!/bin/bash
#SBATCH --job-name=msm8
#SBATCH --partition=fuchs
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=16:00:00
#SBATCH --output=logs/%x_%A_%a.out
# Each array task fits every k-th line of the task file (lines "<CC> <sex>"),
# starting at line SLURM_ARRAY_TASK_ID, where k is the number of array tasks.
# The fuchs QOS allows 20 submitted jobs per user, so the array is shorter than
# the task file. --time covers two countries per task; raise it if there are fewer tasks.
# Submitted by code/hpc_submit.sh, which sets --array and the dependency on compile.sh.
source /home/fuchs/fias/knguyen/.bashrc
conda activate kinh
cd "${SLURM_SUBMIT_DIR}"
TASKS="${1:-code/tasks_f.txt}"
K="${SLURM_ARRAY_TASK_COUNT:-${SLURM_ARRAY_TASK_MAX:-1}}"
I="${SLURM_ARRAY_TASK_ID:-1}"
mapfile -t LINES < <(awk -v k="$K" -v i="$I" 'NF { n++; if ((n - 1) % k == i - 1) print }' "$TASKS")
echo "task ${I} of ${K} on $(hostname): ${#LINES[@]} fits"
for line in "${LINES[@]}"; do
  read -r CC SEX <<< "$line"
  echo "start ${CC} ${SEX} $(date)"
  srun Rscript code/run.r "$CC" "$SEX" < /dev/null || echo "failed ${CC} ${SEX}"
done