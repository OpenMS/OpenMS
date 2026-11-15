#!/bin/bash
# Descriptive diagnosis of the raw02 gap (PLAN_SINGLE_BINARY acceptance verdict unchanged): C++ scorer, seeds 0-4, on
# arm-E arrays with one block replaced by the OpenDIAWorkflow block. Arms: WF (all workflow blocks), F12 (workflow f12),
# D (workflow deficit). Baseline = step 1 (arm-E arrays + featC) in gap-20261114/cppscore.
export LD_LIBRARY_PATH=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/src/tests/class_tests/bin/OpenDIAScorer_arrays
G=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114
cd /home/singjc/.jcode/scratch/e2e/swap
for r in raw02 raw01; do
  A=$G/armE_$r
  for s in 0 1 2 3 4; do
    echo "JCODE_PROGRESS {\"message\":\"$r seed $s\"}"
    $B $A ${r}_WF_s$s.npy $s 8 ${r}_wf_H.npy ${r}_wf_F12.npy ${r}_wf_D.npy ${r}_wf_I.npy > ${r}_WF_s$s.log 2>&1
    $B $A ${r}_F12_s$s.npy $s 8 $A/heuristics.npy ${r}_wf_F12.npy $A/featD.npy $G/cpp/featI_$r.npy > ${r}_F12_s$s.log 2>&1
    $B $A ${r}_D_s$s.npy $s 8 $A/heuristics.npy $A/feat12.npy ${r}_wf_D.npy $G/cpp/featI_$r.npy > ${r}_D_s$s.log 2>&1
  done
done
echo DONE
