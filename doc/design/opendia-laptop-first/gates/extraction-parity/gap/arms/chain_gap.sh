#!/bin/bash
# PLAN_GAP chain, run inside the PAIRLIB SIF with the gap run root bound onto the PAIRLIB path. The run root holds arm E's
# data (lib, lib12, raw0X inputs; hard links, read only) and its own score/ outputs. Usage: chain_gap.sh
set -o pipefail
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
N=$B/pairlib-20261113
S=$N/scripts_gap
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python
ARMS="D11FD D11FDA D11FDG D11FDI D11FDAGI"
WT="D11FDWT D11FDAWT D11FDGWT D11FDIWT D11FDAGIWT"
st() { date -u +"%FT%TZ $1" >> $N/STATUS_GAP.log; }
fail() { st "FAIL $1"; exit 1; }
cd $N || exit 1
sha256sum -c GAP_FROZEN.sha256 > gap.check 2>&1 || fail frozen-hash
mkdir -p score
st start
for r in raw01 raw02; do
  st features-$r
  /usr/bin/time -v $PY $S/features_A.py $r $N > $r/featA.log 2>&1 || fail featA-$r
  /usr/bin/time -v $PY $S/features_G.py $r $N $N/lib12/library12_pop.pqp $N/$r/pipeline/search.features.tsv > $r/featG.log 2>&1 || fail featG-$r
  /usr/bin/time -v $PY $S/features_I.py $r $N $N/lib12/library12_pop.pqp $N/iso_windows.tsv $N/$r/pipeline/search.features.tsv > $r/featI.log 2>&1 || fail featI-$r
  st train-$r
  TRAIN_THREADS=10 OMP_NUM_THREADS=10 /usr/bin/time -v $PY $S/train_gap.py $r $ARMS > score/train_$r.log 2>&1 || fail train-$r
  for a in $ARMS; do
    $PY $N/scripts/make_w.py $r $a lin_$a ${a}W >> score/make_w.log 2>&1 || fail w-$a-$r
    $PY $N/scripts/calib_tail.py $r ${a}W ${a}WT >> score/calib_tail.log 2>&1 || fail t-$a-$r
  done
done
st eval
$PY $N/scripts/evaluate_v2.py $WT > score/EVALUATION.txt 2> score/eval.err || fail eval
$PY $S/compare_gap.py $WT > score/COMPARE.txt 2>> score/eval.err || fail compare
$PY $N/scripts/null_auc.py $WT > score/NULL.txt 2>> score/eval.err || fail null
cd $N/scripts
for a in $WT; do $PY diagnose.py $a $N xic > $N/score/DIAG_$a.txt 2>&1 || fail diag-$a; done
cd $N
for a in $WT; do $PY $S/pooled_rho.py $N $a > score/RHO_$a.txt 2>&1 || fail rho-$a; done
# Control reproduction: this run's D11FDWT against arm E's (deterministic training on the same inputs).
$PY -c "
import numpy as np
for r in ('raw01', 'raw02'):
    a = np.load('$N/' + r + '/scores_D11FDWT.npy'); b = np.load('$N/ref_armE/' + r + '_scores_D11FDWT.npy')
    print(r, 'control max abs diff vs arm E', float(np.abs(a.astype(float) - b.astype(float)).max()))" > score/CONTROL.txt 2>&1 || fail control
st DONE
