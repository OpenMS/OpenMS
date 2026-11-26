#!/usr/bin/env bash
set -euo pipefail
SRC=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment/src/topp/OpenDIA
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${1:-/home/singjc/.jcode/scratch/calibration-selection-replay-20261126}
mkdir -p "$OUT"
for name in selection_isolation selection_real; do
 c++ -O2 -std=c++17 -fopenmp -I"$SRC/include" "$HERE/$name.cpp" "$SRC/source/OpenDIAGBM.cpp" "$SRC/source/OpenDIAScorer.cpp" -o "$OUT/$name"
done
"$OUT/selection_isolation" | tee "$OUT/selection_isolation.txt"
PY=/home/singjc/Documents/github/pyprophet/.venv/bin/python
LIB=/home/singjc/.jcode/scratch/blend-20261125/eval/hela_r5.predicted_decoys.pqp
REF=/home/singjc/.jcode/scratch/blend-20261125/eval/ref/entlib-gbt-20261122/runs
for set in D H; do
 "$PY" "$HERE/export_selection.py" "$LIB" "$REF/${set}_cur_s0/"*.precursors.tsv "$OUT/${set}_net_s0.tsv"
 "$OUT/selection_real" "$OUT/${set}_net_s0.tsv" | tee "$OUT/${set}_selection.txt"
done
"$PY" "$HERE/fold_check.py" | tee "$OUT/fold_check.txt"
