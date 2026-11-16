#!/bin/bash
# Probe (option 1 precheck): C++ CrossFitScorer on the arm-E arrays + C++ featI, with and without 5 signal-level
# interaction inputs (make_interactions.py), seeds 0-2, both raws.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114 || exit 1
export LD_LIBRARY_PATH=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/src/tests/class_tests/bin/OpenDIAScorer_arrays
J=/home/singjc/.jcode/scratch/e2e/joint
mkdir -p $J/score
for r in raw01 raw02; do
  for s in 0 1 2; do
    $B armE_$r $J/score/X_${r}_s$s.npy $s 8 armE_$r/heuristics.npy armE_$r/feat12.npy armE_$r/featD.npy cpp/featI_$r.npy \
      $J/inter_$r.npy > $J/score/X_${r}_s$s.log 2>&1
    tail -1 $J/score/X_${r}_s$s.log
  done
done
