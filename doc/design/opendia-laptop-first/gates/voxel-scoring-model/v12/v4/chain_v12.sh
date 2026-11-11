#!/bin/bash
# PLAN_V12: compute corrected profile views for raw02, then train the dual-view voxel arm on raw01/raw02.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
sha256sum -c v4/PLAN_V12.md.sha256 > v4/PLAN_V12.check 2>&1 || { echo plan-changed > v4/STATUS_V12; exit 1; }
echo features-raw02 > v4/STATUS_V12
[ -f raw02/xic_corrected.npy ] || $PY v4/interference_profiles.py raw02 > v4/interference_raw02.log 2>&1 || { echo failed-features-raw02 > v4/STATUS_V12; exit 1; }
for r in raw01 raw02; do echo train-$r > v4/STATUS_V12
  TRAIN_THREADS=10 OMP_NUM_THREADS=10 /usr/bin/time -v $PY v4/train_v4.py $r D11FI > v4/train_v12_$r.log 2>&1 || { echo failed-train-$r > v4/STATUS_V12; exit 1; }
  $PY v4/make_w.py $r D11FI lin_D11FI D11FIW >> v4/make_w_v12.log 2>&1 || { echo failed-w-$r > v4/STATUS_V12; exit 1; }
  $PY v4/calib_tail.py $r D11FIW D11FIWT >> v4/calib_tail_v12.log 2>&1 || { echo failed-tail-$r > v4/STATUS_V12; exit 1; }
done
echo evaluate > v4/STATUS_V12
$PY evaluate_v2.py D11FWT D11FIWT > v4/EVALUATION_V12.txt 2> v4/eval_v12.err || { echo failed-eval > v4/STATUS_V12; exit 1; }
$PY v3seq/null_auc.py D11FWT D11FIWT > v4/NULL_V12.txt 2>> v4/eval_v12.err
$PY v3seq/misses.py D11FIWT > v4/funnel_D11FIWT.txt 2>> v4/eval_v12.err
$PY v4/vs_base.py D11FWT D11FIWT > v4/VS_D11FWT_V12.txt 2>> v4/eval_v12.err
echo done > v4/STATUS_V12
