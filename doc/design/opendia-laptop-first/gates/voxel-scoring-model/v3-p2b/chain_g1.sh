#!/bin/bash
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
for r in raw01 raw02; do echo g1-$r > v3p2b/STATUS; $PY v3p2b/train_p2b.py $r D11 D11n > v3p2b/train_g1_$r.log 2>&1 || { echo failed-g1-$r > v3p2b/STATUS; exit 1; }; done
echo g1-eval > v3p2b/STATUS
$PY evaluate_v2.py L11 L11d D11 D11n > v3p2b/EVALUATION_G1.txt 2> v3p2b/eval_g1.err || { echo failed-g1-eval > v3p2b/STATUS; exit 1; }
echo g1-done > v3p2b/STATUS
