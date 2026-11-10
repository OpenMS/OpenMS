#!/bin/bash
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
for r in raw01 raw02; do echo g2-$r > v3p2b/STATUS; /usr/bin/time -v $PY v3p2b/train_p2b.py $r DV > v3p2b/train_g2_$r.log 2>&1 || { echo failed-g2-$r > v3p2b/STATUS; exit 1; }; done
echo g2-eval > v3p2b/STATUS
$PY evaluate_v2.py D11 DV FS > v3p2b/EVALUATION_G2.txt 2> v3p2b/eval_g2.err || { echo failed-g2-eval > v3p2b/STATUS; exit 1; }
echo g2-done > v3p2b/STATUS
