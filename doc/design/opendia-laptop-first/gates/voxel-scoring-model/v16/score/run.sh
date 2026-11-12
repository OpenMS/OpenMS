#!/bin/bash
set -uo pipefail
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2/v4/stageD_v16/score
export PYTHONDONTWRITEBYTECODE=1 TRAIN_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python
trap 'rc=$?; echo "$rc" > EXIT; if [ "$rc" = 0 ]; then echo reconstruction-complete > STATUS; else echo failed > STATUS; fi' EXIT
set -e
date -u +%FT%TZ > freeze_reconstruction.timestamp
sha256sum reconstruct.py helpers.py test_helpers.py run.sh ../../train_v4.py ../PLAN_V16.md > freeze_reconstruction.sha256
sha256sum ../../../raw0{1,2}/{scores_D11F.npy,lin_D11F.npy,scores_D11FW.npy,scores_D11FWT.npy,population.tsv,heuristics.npy,feat12.npy,train_windows.tsv} > inputs_before.sha256
$PY test_helpers.py > tests.log 2>&1
for raw in raw01 raw02; do
  echo "reconstructing-$raw" > STATUS
  date -u +%FT%TZ > "freeze_${raw}.timestamp"
  sha256sum reconstruct.py helpers.py ../../train_v4.py > "freeze_${raw}.sha256"
  $PY -u reconstruct.py "$raw" > "$raw.log" 2>&1
done
sha256sum -c inputs_before.sha256 > inputs_after.check
