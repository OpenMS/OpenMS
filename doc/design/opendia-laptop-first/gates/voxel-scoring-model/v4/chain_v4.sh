#!/bin/bash
# PLAN_V4 phase 1: D11R, D11M, D11MR on both raws (10 torch threads), then evaluation, null check, funnel.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
for r in raw01 raw02; do echo train-$r > v4/STATUS
  /usr/bin/time -v $PY v4/train_v4.py $r D11R D11M D11MR > v4/train_$r.log 2>&1 || { echo failed-$r > v4/STATUS; exit 1; }
done
echo evaluate > v4/STATUS
$PY evaluate_v2.py D11 D11R D11M D11MR > v4/EVALUATION_V4.txt 2> v4/eval.err || { echo failed-eval > v4/STATUS; exit 1; }
$PY v3seq/null_auc.py D11 D11R D11M D11MR > v4/NULL_V4.txt 2>> v4/eval.err || { echo failed-null > v4/STATUS; exit 1; }
for a in D11R D11M D11MR; do $PY v3seq/misses.py $a > v4/funnel_$a.txt 2>> v4/eval.err; done
echo done > v4/STATUS
