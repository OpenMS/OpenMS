#!/bin/bash
# PRECHECK_MS1 part 2: C++ CrossFitScorer on the arm-E arrays + C++ featI + 5 MS1 features (ms1_precheck.py), seeds 0-4.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114 || exit 1
export LD_LIBRARY_PATH=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/src/tests/class_tests/bin/OpenDIAScorer_arrays
M=/home/singjc/.jcode/scratch/e2e/ms1
mkdir -p $M/score
for r in raw01 raw02; do
  for s in 0 1 2 3 4; do
    $B armE_$r $M/score/${r}_s$s.npy $s 8 armE_$r/heuristics.npy armE_$r/feat12.npy armE_$r/featD.npy cpp/featI_$r.npy \
      $M/$r/ms1feat.npy > $M/score/${r}_s$s.log 2>&1
    echo "JCODE_PROGRESS {\"message\":\"$r s$s $(tail -1 $M/score/${r}_s$s.log)\"}"
  done
done
