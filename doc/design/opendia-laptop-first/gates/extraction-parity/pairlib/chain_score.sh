#!/bin/bash
# PLAN_PAIRLIB scoring chain (after both chain_raw.sh jobs reach data-done). Same recipe as xicids/interf.
set -o pipefail
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
N=$B/pairlib-20261113
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python
st() { date -u +"%FT%TZ $1" >> $N/STATUS.log; }
fail() { st "FAIL $1"; exit 1; }
cd $N; mkdir -p score
for r in raw01 raw02; do grep -q "$r data-done" STATUS.log || fail "$r-not-done"; done
for r in raw01 raw02; do
  st train-$r
  TRAIN_THREADS=10 OMP_NUM_THREADS=10 /usr/bin/time -v $PY scripts/train_v4.py $r D11F D11FD > score/train_$r.log 2>&1 || fail train-$r
  $PY scripts/make_w.py $r D11F lin_D11F D11FW >> score/make_w.log 2>&1 || fail w-$r
  $PY scripts/calib_tail.py $r D11FW D11FWT >> score/calib_tail.log 2>&1 || fail t-$r
  $PY scripts/make_w.py $r D11FD lin_D11FD D11FDW >> score/make_w.log 2>&1 || fail wd-$r
  $PY scripts/calib_tail.py $r D11FDW D11FDWT >> score/calib_tail.log 2>&1 || fail td-$r
done
st eval
$PY scripts/evaluate_v2.py D11FWT D11FDWT > score/EVALUATION.txt 2> score/eval.err || fail eval
$PY scripts/null_auc.py BASE D11FWT D11FDWT > score/NULL.txt 2>> score/eval.err || fail null
st diag
cd scripts
for a in D11FWT D11FDWT; do
  $PY diagnose.py $a $N xic > $N/score/DIAG_$a.txt 2>&1 || fail diag-$a
  $PY parent.py $a > $N/score/PARENT_$a.txt 2>&1 || fail parent-$a
done
st DONE
