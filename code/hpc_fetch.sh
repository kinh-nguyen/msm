#!/usr/bin/env bash
# Run on your Mac: copies fits, figures and logs back from the cluster.
#   bash code/hpc_fetch.sh
set -euo pipefail
HPC_HOST="${HPC_HOST:-fuchs}"
HPC_DIR="${HPC_DIR:-/panfs/vdura1/fuchs/fias/knguyen/MultistageSurv}"
cd "$(git rev-parse --show-toplevel)"
mkdir -p fit fig logs/hpc
rsync -av "$HPC_HOST:$HPC_DIR/fit/fit8_*" fit/ || true
rsync -av "$HPC_HOST:$HPC_DIR/fig/fit8_*" fig/ || true
rsync -av "$HPC_HOST:$HPC_DIR/logs/" logs/hpc/
grep -h "^done\|Error" logs/hpc/msm8_*.out 2>/dev/null | sort | uniq -c | sort -rn | head -50 || true
