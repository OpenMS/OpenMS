#!/bin/bash
set -euo pipefail
# Execute only after PRECHECK_CALIB commit is pushed.
cd /nfs/research/juan/singj/opendia/calib-20261123
mkdir -p results
PY=("$HOME/.local/bin/uv" run -q --no-project --with pandas --with matplotlib python)
for P in D H; do
  "${PY[@]}" eval_calib.py /nfs/research/juan/singj/opendia/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp /nfs/research/juan/singj/opendia/gbmnull-20261123/lay/hela "$P" "results/$P"
done
for P in A1 A2 A3; do
  "${PY[@]}" eval_calib.py /nfs/research/juan/singj/opendia/gap-20261114/home/lib/library.pqp /nfs/research/juan/singj/opendia/gbmnull-20261123/lay/astral "$P" "results/$P"
done
"${PY[@]}" plot_calib.py results
