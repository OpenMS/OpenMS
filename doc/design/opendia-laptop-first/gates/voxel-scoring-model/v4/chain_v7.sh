#!/bin/bash
# PLAN_V7: D11FX trained (train_v4.py), W rule -> D11FXW, tail transform -> D11FXWT; evaluation against BASE and D11FWT.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
for r in raw01 raw02; do echo train-$r > v4/STATUS_V7
  TRAIN_THREADS=4 OMP_NUM_THREADS=4 /usr/bin/time -v $PY v4/train_v4.py $r D11FX > v4/train_v7_$r.log 2>&1 || { echo failed-$r > v4/STATUS_V7; exit 1; }
  $PY v4/make_w.py $r D11FX lin_D11FX D11FXW >> v4/make_w_v7.log 2>&1 || { echo failed-w > v4/STATUS_V7; exit 1; }
  $PY v4/calib_tail.py $r D11FXW D11FXWT >> v4/calib_tail_v7.log 2>&1 || { echo failed-t > v4/STATUS_V7; exit 1; }
done
echo evaluate > v4/STATUS_V7
$PY evaluate_v2.py D11FWT D11FXWT > v4/EVALUATION_V7.txt 2> v4/eval_v7.err || { echo failed-eval > v4/STATUS_V7; exit 1; }
$PY v3seq/null_auc.py D11FWT D11FXWT > v4/NULL_V7.txt 2>> v4/eval_v7.err
$PY v3seq/misses.py D11FXWT > v4/funnel_D11FXWT.txt 2>> v4/eval_v7.err
$PY v4/vs_base.py D11FWT D11FXWT > v4/VS_D11FWT_V7.txt 2>> v4/eval_v7.err
echo done > v4/STATUS_V7
