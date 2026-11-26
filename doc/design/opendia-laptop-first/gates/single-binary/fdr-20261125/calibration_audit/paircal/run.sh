#!/usr/bin/env bash
set -euo pipefail
s=/home/singjc/.jcode/scratch/paircal-20261126
p=/home/singjc/.jcode/scratch/mlenv/bin/python
lib=/home/singjc/.jcode/scratch/blend-20261125/eval/hela_r5.predicted_decoys.pqp
"$p" "$s/paircal_helpers.py"
"$p" -u "$s/run_selection.py" > "$s/selection.log" 2>&1
gzip -1f "$s/selected.tsv"
"$p" -u "$s/evaluate.py" "$lib" "$s/selected.tsv.gz" "$s/evaluation.tsv" > "$s/evaluation.txt" 2>&1
