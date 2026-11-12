#!/bin/bash
# PLAN_V18: featS -> control + D11FS training -> W -> T -> evaluation. Run from v2/v4/v18.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2/v4/v18
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python; V2=../..
trap 'rc=$?; echo $rc > EXIT' EXIT
sha256sum -c PLAN_V18.sha256 > PLAN_V18.check 2>&1 || { echo plan-changed > STATUS; exit 1; }
$PY test_sib_core.py > tests.log 2>&1 || { echo failed-tests > STATUS; exit 1; }
for r in raw01 raw02; do
  echo features-$r > STATUS; $PY build_sib.py $r >> features.log 2>&1 || { echo failed-features > STATUS; exit 1; }
  echo train-$r > STATUS
  OMP_NUM_THREADS=10 /usr/bin/time -v $PY -u run_v18.py $r > train_$r.log 2>&1 || { echo failed-train-$r > STATUS; exit 1; }
  $PY ../make_w.py $r D11FS lin_D11FS D11FSW >> make_w.log 2>&1 || { echo failed-w > STATUS; exit 1; }
  $PY ../calib_tail.py $r D11FSW D11FSWT >> calib_tail.log 2>&1 || { echo failed-t > STATUS; exit 1; }
done
echo evaluate > STATUS
cd $V2
$PY evaluate_v2.py D11FWT D11FSWT > v4/v18/EVALUATION_V18.txt 2> v4/v18/eval.err || { echo failed-eval > v4/v18/STATUS; exit 1; }
$PY v3seq/null_auc.py D11FWT D11FSWT > v4/v18/NULL_V18.txt 2>> v4/v18/eval.err
$PY v4/vs_base.py D11FWT D11FSWT > v4/v18/VS_D11FWT_V18.txt 2>> v4/v18/eval.err
$PY v3seq/misses.py D11FSWT > v4/v18/funnel_D11FSWT.txt 2>> v4/v18/eval.err
cd v4/v18 && $PY gains.py D11FWT D11FSWT > GAINS_V18.txt 2>> eval.err
echo done > STATUS
