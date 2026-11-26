#!/usr/bin/env bash
# Pooled arm A and AB (replay) evaluation on finished A dev runs, side by side with G0 (= B0 replay = exported tool score).
set -uo pipefail
H=/home/singjc/.jcode/scratch/fdr-ab
LIB=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114/armE_raw01/library.pqp
PY=/home/singjc/Documents/github/pyprophet/.venv/bin/python
mkdir -p $H/screenAB
A=(); AB=()
for s in 0 1 2 3 4; do for r in 01 02; do
  f=$H/runs/A_raw${r}_s$s/LFQ_Astral_DIA_15min_50ng_Human_${r}.features.tsv
  [ -s $H/runs/A_raw${r}_s$s/precursors.tsv ] || continue
  A+=("$f:self")
  o=$H/screenAB/AB2_raw${r}_s$s.tsv
  [ -s $o ] || nice $H/gbm_replay $f $H/astral_keys.tsv $s 6 $o mono=i_explained:-1,i_nclaim:-1,e_explained:-1,e_nclaim:-1 > ${o%.tsv}.log 2>&1
  AB+=("$f:$o")
done; done
$PY $H/eval_replay.py $LIB A "${A[@]}" > $H/screenAB/EVAL_A.txt
$PY $H/eval_replay.py $LIB AB2 "${AB[@]}" > $H/screenAB/EVAL_AB2.txt
cat $H/screenAB/EVAL_A.txt $H/screenAB/EVAL_AB2.txt
