#!/bin/bash
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
for r in raw01 raw02; do echo train-$r > v3p1/STATUS_P2a; $PY v3p1/train_p1.py $r FL11 G11a > v3p1/train_p2a_$r.log 2>&1 || { echo failed-$r > v3p1/STATUS_P2a; exit 1; }; done
echo evaluate > v3p1/STATUS_P2a
$PY evaluate_v2.py L11 G11 G11a FL11 > v3p1/EVALUATION_P2a.txt 2> v3p1/evaluate_p2a.err || { echo eval-failed > v3p1/STATUS_P2a; exit 1; }
echo done > v3p1/STATUS_P2a
