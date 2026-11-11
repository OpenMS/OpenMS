#!/bin/bash
# PLAN_V5: D11FP, D11FPC trained (train_v4.py), W rule -> D11FPW, D11FPCW; evaluation against BASE, D11FW, D11.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
for r in raw01 raw02; do echo train-$r > v4/STATUS_V5
  OMP_NUM_THREADS=10 /usr/bin/time -v $PY v4/train_v4.py $r D11FP D11FPC > v4/train_v5_$r.log 2>&1 || { echo failed-$r > v4/STATUS_V5; exit 1; }
  $PY v4/make_w.py $r D11FP lin_D11FP D11FPW >> v4/make_w_v5.log 2>&1 && $PY v4/make_w.py $r D11FPC lin_D11FPC D11FPCW >> v4/make_w_v5.log 2>&1 || { echo failed-w > v4/STATUS_V5; exit 1; }
done
echo evaluate > v4/STATUS_V5
$PY evaluate_v2.py D11 D11FW D11FPW D11FPCW > v4/EVALUATION_V5.txt 2> v4/eval_v5.err || { echo failed-eval > v4/STATUS_V5; exit 1; }
$PY v3seq/null_auc.py D11FW D11FPW D11FPCW > v4/NULL_V5.txt 2>> v4/eval_v5.err
for a in D11FPW D11FPCW; do $PY v3seq/misses.py $a > v4/funnel_$a.txt 2>> v4/eval_v5.err; done
$PY v4/vs_base.py D11FW D11FPW D11FPCW > v4/VS_D11FW_V5.txt 2>> v4/eval_v5.err
$PY v4/vs_base.py D11 D11FPW D11FPCW > v4/VS_D11_V5.txt 2>> v4/eval_v5.err
echo done > v4/STATUS_V5
