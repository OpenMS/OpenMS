#!/usr/bin/env bash
# Local evaluation only, after fetching the complete remote rtnull-20261126 workspace.
# Usage: bash rtnull_eval.sh /absolute/path/to/fetched/rtnull-20261126
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
FETCH=$1
PY=/home/singjc/Documents/github/pyprophet/.venv/bin/python
LIB=/home/singjc/.jcode/scratch/blend-20261125/eval/hela_r5.predicted_decoys.pqp
REF=/home/singjc/.jcode/scratch/blend-20261125/eval/ref/entlib-gbt-20261122/runs
LAY=$FETCH/eval_layout
mkdir "$LAY"
failed=0
for P in D H; do
 for S in 0 1 2 3 4; do
  O=$FETCH/runs/${P}_NORT_s$S
  if [[ ! -f "$O.exit" ]] || [[ $(cat "$O.exit") != 0 ]] || [[ $(find "$O" -maxdepth 1 -name '*.precursors.tsv' 2>/dev/null | wc -l) != 1 ]]; then
   echo "INCOMPLETE_OR_FAILED $P seed $S, retained in fixed denominator"
   failed=1
  fi
  ln -s "$REF/${P}_gbt_s$S" "$LAY/${P}_G0_s$S"
  ln -s "$O" "$LAY/${P}_NORT_s$S"
 done
done
if [[ $failed != 0 ]]; then echo 'No pass evaluable: not all ten frozen cells succeeded.'; exit 3; fi
sha256sum -c "$HERE/baseline.sha256"
for P in D H; do
 RELABEL=$HERE/label_audit_hela.tsv "$PY" "$HERE/eval_pdec_relabel.py" "$LIB" "$LAY" "$P"
done
