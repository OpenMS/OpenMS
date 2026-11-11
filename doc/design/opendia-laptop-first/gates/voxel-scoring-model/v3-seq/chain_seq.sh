#!/bin/bash
# PLAN_V3_SEQ: DS11 on both raws, then DVS on both raws (one process at a time, 10 threads); then evaluation, null check, vs control.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python; export OMP_NUM_THREADS=10
for arm in DS11 DVS; do for r in raw01 raw02; do
  echo $arm-$r > v3seq/STATUS
  /usr/bin/time -v $PY v3seq/train_seq.py $r $arm > v3seq/train_${arm}_$r.log 2>&1 || { echo failed-$arm-$r > v3seq/STATUS; exit 1; }
done; done
echo evaluate > v3seq/STATUS
$PY evaluate_v2.py D11 DS11 DV DVS > v3seq/EVALUATION_SEQ.txt 2> v3seq/eval.err || { echo failed-eval > v3seq/STATUS; exit 1; }
$PY v3seq/null_auc.py D11 DS11 DV DVS > v3seq/NULL_SEQ.txt 2>> v3seq/eval.err || { echo failed-null > v3seq/STATUS; exit 1; }
$PY v3seq/vs_control.py > v3seq/VS_CONTROL.txt 2>> v3seq/eval.err || { echo failed-vs > v3seq/STATUS; exit 1; }
echo done > v3seq/STATUS
