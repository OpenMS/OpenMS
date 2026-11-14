#!/bin/bash
# PLAN_SINGLE_BINARY step 1: C++ CrossFitScorer over seeds 0-4 on the arm-E arrays + C++ featI, both raws.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114 || exit 1
export LD_LIBRARY_PATH=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/src/tests/class_tests/bin/OpenDIAScorer_arrays
for r in raw01 raw02; do
  for s in 0 1 2 3 4; do
    $B armE_$r cppscore/${r}_s$s.npy $s 8 armE_$r/heuristics.npy armE_$r/feat12.npy armE_$r/featD.npy cpp/featI_$r.npy \
      > cppscore/${r}_s$s.log 2>&1
  done
done
