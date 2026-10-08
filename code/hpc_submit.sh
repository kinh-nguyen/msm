#!/usr/bin/env bash
# Run on your Mac, from anywhere inside the project:
#   bash code/hpc_submit.sh                 all 37 countries, women
#   bash code/hpc_submit.sh code/tasks_m.txt   men
#   HPC_HOST=user@host bash code/hpc_submit.sh  if "fuchs" is not an ssh alias
# Syncs the pipeline and data to the cluster, compiles once, then submits one
# array task per line of the task file. Needs data/cc8.csv.bz2 (rerun prep.Rmd first).
set -euo pipefail
HPC_HOST="${HPC_HOST:-fuchs}"
HPC_DIR="${HPC_DIR:-/panfs/vdura1/fuchs/fias/knguyen/MultistageSurv}"
TASKS="${1:-code/tasks_f.txt}"
MAXRUN="${MAXRUN:-10}" # array tasks running at once

cd "$(git rev-parse --show-toplevel)"
[ -f data/cc8.csv.bz2 ] || { echo "data/cc8.csv.bz2 missing: rerun code/prep.Rmd first"; exit 1; }
N=$(grep -c . "$TASKS")

ssh "$HPC_HOST" "mkdir -p '$HPC_DIR/code' '$HPC_DIR/data' '$HPC_DIR/fit' '$HPC_DIR/fig' '$HPC_DIR/logs' && touch '$HPC_DIR/.here'"
rsync -av code/msm.cpp code/ktools.hpp code/run.r code/compile.r code/run.sh code/compile.sh "$TASKS" \
  "$HPC_HOST:$HPC_DIR/code/"
rsync -av data/cc8.csv.bz2 "$HPC_HOST:$HPC_DIR/data/"
for f in data/pd.csv.bz2 data/mjs.csv; do [ -f "$f" ] && rsync -av "$f" "$HPC_HOST:$HPC_DIR/data/"; done

ssh "$HPC_HOST" bash -s <<REMOTE
set -e
cd '$HPC_DIR'
jid=\$(sbatch --parsable code/compile.sh)
aid=\$(sbatch --parsable --dependency=afterok:\$jid --array=1-$N%$MAXRUN code/run.sh code/$(basename "$TASKS"))
echo "compile job \$jid, array job \$aid ($N tasks)"
REMOTE
echo 'check progress with: ssh '"$HPC_HOST"' squeue -u $USER'
