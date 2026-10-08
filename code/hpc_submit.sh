#!/usr/bin/env bash
# Run on your Mac, from anywhere inside the project:
#   bash code/hpc_submit.sh                 all 37 countries, women
#   bash code/hpc_submit.sh code/tasks_m.txt   men
#   HPC_HOST=user@host bash code/hpc_submit.sh  if "fuchs" is not an ssh alias
# Syncs the pipeline and data to the cluster, compiles once, then submits an
# array whose tasks share the lines of the task file. Needs data/cc8.csv.bz2.
# The fuchs QOS allows MAXSUB = 20 submitted and 10 running jobs per user; the
# array length is reduced to fit beside the compile job and any jobs already queued.
set -euo pipefail
HPC_HOST="${HPC_HOST:-knguyen@fuchs.hhlr-gu.de}"
HPC_DIR="${HPC_DIR:-/panfs/vdura1/fuchs/fias/knguyen/MultistageSurv}"
TASKS="${1:-code/tasks_f.txt}"
MAXRUN="${MAXRUN:-10}" # array tasks running at once (QOS MaxJobsPU)
MAXSUB="${MAXSUB:-20}" # jobs submitted at once (QOS MaxSubmitPU)

cd "$(git rev-parse --show-toplevel)"
[ -f data/cc8.csv.bz2 ] || { echo "data/cc8.csv.bz2 missing: rerun code/prep.Rmd first"; exit 1; }
N=$(grep -c . "$TASKS")

ssh "$HPC_HOST" "mkdir -p '$HPC_DIR/code' '$HPC_DIR/data' '$HPC_DIR/fit' '$HPC_DIR/fig' '$HPC_DIR/logs' && touch '$HPC_DIR/.here'"
rsync -av code/msm.cpp code/run.r code/compile.r code/run.sh code/compile.sh "$TASKS" \
  "$HPC_HOST:$HPC_DIR/code/"
rsync -av data/cc8.csv.bz2 "$HPC_HOST:$HPC_DIR/data/"
for f in data/pd.csv.bz2 data/mjs.csv; do [ -f "$f" ] && rsync -av "$f" "$HPC_HOST:$HPC_DIR/data/"; done

ssh "$HPC_HOST" bash -s <<REMOTE
set -e
cd '$HPC_DIR'
queued=\$(squeue -u \$USER -h -r | wc -l)
free=\$(( $MAXSUB - queued - 1 ))
[ \$free -ge 1 ] || { echo "\$queued jobs already queued, QOS allows $MAXSUB"; exit 1; }
narray=\$(( $N < free ? $N : free ))
jid=\$(sbatch --parsable code/compile.sh)
aid=\$(sbatch --parsable --dependency=afterok:\$jid --array=1-\$narray%$MAXRUN code/run.sh code/$(basename "$TASKS"))
echo "compile job \$jid, array job \$aid (\$narray tasks for $N fits)"
REMOTE
echo 'check progress with: ssh '"$HPC_HOST"' squeue -u $USER'