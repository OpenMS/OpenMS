#!/bin/bash
# PLAN_V9: features N for raw 02 (raw 01 already built), D11FN trained, W rule -> D11FNW, tail transform -> D11FNWT;
# evaluation against BASE and D11FWT, null AUC, funnel, paired bootstrap vs D11FWT.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
sha256sum -c v4/PLAN_V9.md.sha256 > v4/PLAN_V9.check 2>&1 || { echo plan-changed > v4/STATUS_V9; exit 1; }
echo features-raw02 > v4/STATUS_V9
[ -f raw02/nat.npy ] || $PY v4/native_features.py raw02 > v4/native_raw02.log 2>&1 || { echo failed-features > v4/STATUS_V9; exit 1; }
for r in raw01 raw02; do echo train-$r > v4/STATUS_V9
  TRAIN_THREADS=10 OMP_NUM_THREADS=10 /usr/bin/time -v $PY v4/train_v4.py $r D11FN > v4/train_v9_$r.log 2>&1 || { echo failed-$r > v4/STATUS_V9; exit 1; }
  $PY v4/make_w.py $r D11FN lin_D11FN D11FNW >> v4/make_w_v9.log 2>&1 || { echo failed-w > v4/STATUS_V9; exit 1; }
  $PY v4/calib_tail.py $r D11FNW D11FNWT >> v4/calib_tail_v9.log 2>&1 || { echo failed-t > v4/STATUS_V9; exit 1; }
done
echo evaluate > v4/STATUS_V9
$PY evaluate_v2.py D11FWT D11FNWT > v4/EVALUATION_V9.txt 2> v4/eval_v9.err || { echo failed-eval > v4/STATUS_V9; exit 1; }
$PY v3seq/null_auc.py D11FWT D11FNWT > v4/NULL_V9.txt 2>> v4/eval_v9.err
$PY v3seq/misses.py D11FNWT > v4/funnel_D11FNWT.txt 2>> v4/eval_v9.err
$PY v4/vs_base.py D11FWT D11FNWT > v4/VS_D11FWT_V9.txt 2>> v4/eval_v9.err
echo done > v4/STATUS_V9
