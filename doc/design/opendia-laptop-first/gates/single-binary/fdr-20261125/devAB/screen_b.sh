#!/usr/bin/env bash
# Offline arm-B screen on Astral dev (raw01/02 x s0-4): replay the tool's GBM with each pre-declared variant.
set -uo pipefail
H=/home/singjc/.jcode/scratch/fdr-ab; D=/home/singjc/.jcode/scratch/fdr-20261125
LIB=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114/armE_raw01/library.pqp
PY=/home/singjc/Documents/github/pyprophet/.venv/bin/python
I=i_explained,i_libcos_resid,i_dlibcos,i_nclaim
declare -A V=(
  [B0]=""
  [B1]="exclude=$I"
  [B2]="mono=i_explained:-1,i_nclaim:-1"
  [B3]="exclude=i_libcos_resid,i_dlibcos mono=i_explained:-1,i_nclaim:-1"
  [B4]="exclude=$I,f12_extra_support,f12_extra_coelution,f12_extra_share"
)
mkdir -p $H/screen
for v in B0 B1 B2 B3 B4; do
  pairs=()
  for r in raw01 raw02; do for s in 0 1 2 3 4; do
    f=$(ls $D/${r}_s$s/*.features.tsv); o=$H/screen/${v}_${r}_s$s.tsv
    [ -s $o ] || $H/gbm_replay $f $H/astral_keys.tsv $s 6 $o ${V[$v]} > $H/screen/${v}_${r}_s$s.log 2>&1 || echo "FAIL $v $r $s"
    pairs+=("$f:$o")
  done; done
  $PY $H/eval_replay.py $LIB $v "${pairs[@]}" | tee -a $H/screen/EVAL.txt
done
