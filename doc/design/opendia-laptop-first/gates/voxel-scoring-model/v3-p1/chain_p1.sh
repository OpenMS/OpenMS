#!/bin/bash
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
for r in raw01 raw02; do echo train-$r > v3p1/STATUS; $PY v3p1/train_p1.py $r L11 N11 N11a > v3p1/train_$r.log 2>&1 || { echo failed-$r > v3p1/STATUS; exit 1; }; done
echo evaluate > v3p1/STATUS
$PY evaluate_v2.py L11 N11 N11a F0 > v3p1/EVALUATION_P1.txt 2> v3p1/evaluate.err || { echo eval-failed > v3p1/STATUS; exit 1; }
echo done > v3p1/STATUS
