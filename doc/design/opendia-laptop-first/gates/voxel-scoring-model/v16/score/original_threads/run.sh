#!/bin/bash
set -euo pipefail
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2/v4/stageD_v16/score/original_threads
export PYTHONDONTWRITEBYTECODE=1 TRAIN_THREADS=10
unset MKL_NUM_THREADS
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python
trap 'rc=$?; echo "$rc" > EXIT; if [ "$rc" = 0 ]; then echo reconstruction-complete > STATUS; else echo failed > STATUS; fi' EXIT
date -u +%FT%TZ > freeze_runner.timestamp
sha256sum reconstruct.py helpers.py test_helpers.py run.sh ../../../train_v4.py ../../PLAN_V16.md > freeze_runner.sha256
sha256sum ../../../../raw0{1,2}/{scores_D11F.npy,lin_D11F.npy,scores_D11FW.npy,scores_D11FWT.npy,population.tsv,heuristics.npy,feat12.npy,train_windows.tsv} > inputs_before.sha256
$PY test_helpers.py > tests.log 2>&1
for raw in raw01 raw02; do
  export OMP_NUM_THREADS=10
  if [ "$raw" = raw01 ]; then export OMP_NUM_THREADS=6; fi
  echo "reconstructing-$raw" > STATUS
  date -u +%FT%TZ > "freeze_${raw}.timestamp"
  sha256sum reconstruct.py helpers.py ../../../train_v4.py > "freeze_${raw}.sha256"
  env | sort | rg '^(OMP|MKL|TRAIN|PYTHONDONT)' > "settings_${raw}.txt"
  $PY -u reconstruct.py "$raw" > "$raw.log" 2>&1
done
sha256sum -c inputs_before.sha256 > inputs_after.check
