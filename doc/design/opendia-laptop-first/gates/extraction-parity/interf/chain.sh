#!/bin/bash
# PLAN_INTERF chain. One run.
set -u
N=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/interf-20261113
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python
cd $N
st() { date -u +"%FT%TZ $1" >> STATUS.log; }
fail() { st "FAIL $1"; exit 1; }
sha256sum -c PLAN_INTERF.md.sha256 scripts.sha256 > plan.check 2>&1 || fail hash
for r in raw01 raw02; do
  st feat-$r; $PY scripts/deficit_features.py $r $N >> score/features.log 2>&1 || fail feat-$r
done
for r in raw01 raw02; do
  st train-$r; $PY scripts/train_v4.py $r D11FD D11F > score/train_$r.log 2>&1 || fail train-$r
  $PY scripts/make_w.py $r D11FD lin_D11FD D11FDW >> score/make_w.log 2>&1 || fail w-$r
  $PY scripts/calib_tail.py $r D11FDW D11FDWT >> score/calib_tail.log 2>&1 || fail t-$r
  $PY scripts/make_w.py $r D11F lin_D11F D11FW_rerun >> score/make_w.log 2>&1 || fail wc-$r
  $PY scripts/calib_tail.py $r D11FW_rerun D11FWT_rerun >> score/calib_tail.log 2>&1 || fail tc-$r
done
st eval; $PY scripts/evaluate_v2.py D11FWT D11FWT_rerun D11FDWT > score/EVALUATION.txt 2> score/eval.err || fail eval
st null; $PY scripts/null_auc.py D11FWT D11FWT_rerun D11FDWT > score/NULL.txt 2>&1 || fail null
st diag; (cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/qcal-20261113 && $PY diagnose.py D11FDWT $N xic > $N/score/DIAG_D11FDWT.txt 2>&1) || fail diag
st DONE
