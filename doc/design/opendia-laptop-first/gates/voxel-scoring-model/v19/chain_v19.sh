#!/bin/bash
# PLAN_V19 chain. Run from v2/v4/v19. raw01 peak.npy was built before the plan (look); raw02 is built here.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2/v4/v19
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python; V2=../..
trap 'rc=$?; echo $rc > EXIT' EXIT
sha256sum -c PLAN_V19.sha256 > PLAN_V19.check 2>&1 || { echo plan-changed > STATUS; exit 1; }
$PY test_peak_core.py > tests.log 2>&1 || { echo failed-tests > STATUS; exit 1; }
echo features-raw02 > STATUS
/usr/bin/time -v $PY -u build_peak.py raw02 > build_raw02.log 2>&1 || { echo failed-features > STATUS; exit 1; }
for r in raw01 raw02; do
  echo train-$r > STATUS
  OMP_NUM_THREADS=10 /usr/bin/time -v $PY -u run_v19.py $r > train_$r.log 2>&1 || { echo failed-train-$r > STATUS; exit 1; }
  $PY ../make_w.py $r D11FG lin_D11FG D11FGW >> make_w.log 2>&1 || { echo failed-w > STATUS; exit 1; }
  $PY ../calib_tail.py $r D11FGW D11FGWT >> calib_tail.log 2>&1 || { echo failed-t > STATUS; exit 1; }
done
echo evaluate > STATUS
cd $V2
$PY evaluate_v2.py D11FWT D11FGWT > v4/v19/EVALUATION_V19.txt 2> v4/v19/eval.err || { echo failed-eval > v4/v19/STATUS; exit 1; }
$PY v3seq/null_auc.py D11FWT D11FGWT > v4/v19/NULL_V19.txt 2>> v4/v19/eval.err
$PY v4/vs_base.py D11FWT D11FGWT > v4/v19/VS_D11FWT_V19.txt 2>> v4/v19/eval.err
$PY v3seq/misses.py D11FGWT > v4/v19/funnel_D11FGWT.txt 2>> v4/v19/eval.err
cd v4/v19 && $PY gains.py D11FWT D11FGWT > GAINS_V19.txt 2>> eval.err
echo done > STATUS
