#!/bin/bash
set -euo pipefail
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2/v4/stageD_v16/score/original_threads
export PYTHONDONTWRITEBYTECODE=1
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python
trap 'rc=$?; echo "$rc" > ANALYSIS_EXIT; if [ "$rc" = 0 ]; then echo complete > ANALYSIS_STATUS; else echo failed > ANALYSIS_STATUS; fi' EXIT
echo waiting-reconstruction > ANALYSIS_STATUS
for i in $(seq 1 240); do
  if [ -f EXIT ]; then break; fi
  sleep 2
done
test -f EXIT
test "$(cat EXIT)" = 0
sha256sum raw0{1,2}/{checkpoint_fold*_r*.pt,membership_fold*_r*.npz,selectors_fold*_r*.npz,LABEL_BLIND_COMPLETE,verification.json} > label_blind_artifacts.sha256
date -u +%FT%TZ > freeze_coverage.timestamp
sha256sum coverage.py mismatch_helpers.py helpers.py test_mismatch.py run_analysis.sh > freeze_coverage.sha256
$PY test_mismatch.py > tests_mismatch.log 2>&1
echo analyzing > ANALYSIS_STATUS
$PY -u coverage.py > coverage.log 2>&1
sha256sum -c inputs_before.sha256 > inputs_final.check
