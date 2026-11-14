#!/bin/bash
# PLAN_GAP amendment 1: the evaluation steps of chain_gap.sh after `compare` (which failed on a variable overwrite in
# compare_gap.py), unchanged otherwise. Scores and EVALUATION.txt are not recomputed.
set -o pipefail
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
N=$B/pairlib-20261113
S=$N/scripts_gap
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python
WT="D11FDWT D11FDAWT D11FDGWT D11FDIWT D11FDAGIWT"
st() { date -u +"%FT%TZ $1" >> $N/STATUS_GAP.log; }
fail() { st "FAIL $1"; exit 1; }
cd $N || exit 1
sha256sum -c GAP_AMEND1.sha256 > gap_amend1.check 2>&1 || fail amend1-hash
st eval-amend1
$PY $S/compare_gap.py $WT > score/COMPARE.txt 2>> score/eval.err || fail compare
$PY $N/scripts/null_auc.py $WT > score/NULL.txt 2>> score/eval.err || fail null
cd $N/scripts
for a in $WT; do $PY diagnose.py $a $N xic > $N/score/DIAG_$a.txt 2>&1 || fail diag-$a; done
cd $N
for a in $WT; do $PY $S/pooled_rho.py $N $a > score/RHO_$a.txt 2>&1 || fail rho-$a; done
$PY -c "
import numpy as np
for r in ('raw01', 'raw02'):
    a = np.load('$N/' + r + '/scores_D11FDWT.npy'); b = np.load('$N/ref_armE/' + r + '_scores_D11FDWT.npy')
    print(r, 'control max abs diff vs arm E', float(np.abs(a.astype(float) - b.astype(float)).max()))" > score/CONTROL.txt 2>&1 || fail control
st DONE
