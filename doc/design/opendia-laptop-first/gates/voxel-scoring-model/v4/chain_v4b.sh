#!/bin/bash
# PLAN_V4B: D11W now (raw01, raw02), then D11F/D11FW per raw once its feat12.npy exists, then evaluation.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python; L=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-lib12-20261111
for r in raw01 raw02; do $PY v4/make_w.py $r D11 scores_D11lin D11W >> v4/make_w.log 2>&1 || { echo failed-w > v4/STATUS_B; exit 1; }; done
# raw02 12-fragment arrays: wait for its export.
while [ "$(cat $L/EXPORT_STATUS)" != "done" ]; do [ "$(cat $L/EXPORT_STATUS | cut -c1-6)" = failed ] && { echo failed-export > v4/STATUS_B; exit 1; }; echo wait-export02 > v4/STATUS_B; sleep 30; done
echo train-raw01 > v4/STATUS_B
OMP_NUM_THREADS=6 /usr/bin/time -v $PY v4/train_v4.py raw01 D11F > v4/train_b_raw01.log 2>&1 &
TP=$!
echo build12-raw02 > v4/STATUS_B
/usr/bin/time -v $PY v4/build12.py raw02 > v4/build12_raw02.log 2>&1 && $PY v4/features12.py raw02 > v4/features12_raw02.log 2>&1 || { echo failed-build02 > v4/STATUS_B; exit 1; }
wait $TP || { echo failed-train01 > v4/STATUS_B; exit 1; }
echo train-raw02 > v4/STATUS_B
OMP_NUM_THREADS=10 /usr/bin/time -v $PY v4/train_v4.py raw02 D11F > v4/train_b_raw02.log 2>&1 || { echo failed-train02 > v4/STATUS_B; exit 1; }
for r in raw01 raw02; do $PY v4/make_w.py $r D11F lin_D11F D11FW >> v4/make_w.log 2>&1; done
echo evaluate > v4/STATUS_B
$PY evaluate_v2.py D11 D11W D11F D11FW > v4/EVALUATION_V4B.txt 2> v4/eval_b.err || { echo failed-eval > v4/STATUS_B; exit 1; }
$PY v3seq/null_auc.py D11 D11W D11F D11FW > v4/NULL_V4B.txt 2>> v4/eval_b.err
for a in D11W D11F D11FW; do $PY v3seq/misses.py $a > v4/funnel_$a.txt 2>> v4/eval_b.err; done
$PY v4/vs_d11.py D11W D11F D11FW > v4/VS_D11_B.txt 2>> v4/eval_b.err
echo done > v4/STATUS_B
