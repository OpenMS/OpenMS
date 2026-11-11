#!/bin/bash
# PLAN_V8: D11F -> W rule -> tail transform (D11FWT) on the rtlin population; evaluation against the rebuilt BASE.
# Waits for chain_data.sh to reach data-done (polls STATUS every 60 s, stops on failed-*).
N=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/rtlin-20261111
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
cd $N
while true; do s=$(cat STATUS); [ "$s" = data-done ] && break; case $s in failed-*) echo "data-$s" > v4/STATUS_V8; exit 1;; esac; sleep 60; done
sha256sum -c v4/PLAN_V8.md.sha256 > v4/PLAN_V8.check 2>&1 || { echo plan-changed > v4/STATUS_V8; exit 1; }
for r in raw01 raw02; do echo train-$r > v4/STATUS_V8
  TRAIN_THREADS=10 OMP_NUM_THREADS=10 /usr/bin/time -v $PY scripts/train_v4.py $r D11F > v4/train_v8_$r.log 2>&1 || { echo failed-$r > v4/STATUS_V8; exit 1; }
  $PY scripts/make_w.py $r D11F lin_D11F D11FW >> v4/make_w_v8.log 2>&1 || { echo failed-w > v4/STATUS_V8; exit 1; }
  $PY scripts/calib_tail.py $r D11FW D11FWT >> v4/calib_tail_v8.log 2>&1 || { echo failed-t > v4/STATUS_V8; exit 1; }
done
echo evaluate > v4/STATUS_V8
$PY scripts/evaluate_v2.py D11FW D11FWT > v4/EVALUATION_V8.txt 2> v4/eval_v8.err || { echo failed-eval > v4/STATUS_V8; exit 1; }
$PY scripts/null_auc.py BASE D11FWT > v4/NULL_V8.txt 2>> v4/eval_v8.err
$PY scripts/misses.py D11FWT > v4/funnel_D11FWT.txt 2>> v4/eval_v8.err
echo done > v4/STATUS_V8
