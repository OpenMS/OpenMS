#!/bin/bash
# Voxel pilot chain: waits for export01 (itself queued behind round 7), checks the export kept the base candidate rows,
# then builds voxels, trains the 8 cells and evaluates under PLAN_PILOT.md. STATUS_PILOT shows the step.
O=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
cd $O
echo waiting-export > STATUS_PILOT
while [ "$(cat export01/STATUS 2>/dev/null)" != "done" ]; do sleep 20; done
if [ "$(cat export01/EXIT)" != "0" ] || ! grep -q IDENTICAL export01/CMP.txt; then echo export-failed > STATUS_PILOT; exit 1; fi
echo voxels > STATUS_PILOT
/usr/bin/time -v $PY build_voxels.py > voxels.log 2> voxels.err || { echo voxels-failed > STATUS_PILOT; exit 1; }
echo train > STATUS_PILOT
/usr/bin/time -v $PY train_pilot.py > train.log 2> train.err || { echo train-failed > STATUS_PILOT; exit 1; }
echo evaluate > STATUS_PILOT
$PY evaluate_pilot.py > EVALUATION_PILOT.txt 2> evaluate.err || { echo evaluate-failed > STATUS_PILOT; exit 1; }
echo done > STATUS_PILOT
