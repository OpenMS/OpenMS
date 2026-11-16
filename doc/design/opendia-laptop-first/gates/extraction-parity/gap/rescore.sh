#!/bin/bash
# PRECHECK_PEAK: rerun the C++ scorer (seeds 0-2) on the arm-E arrays + C++ featI, writing every window's raw score and
# linear part, so the window-choice what-ifs can be evaluated offline. Scores equal cppscore/ (same code, same seeds).
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114 || exit 1
export LD_LIBRARY_PATH=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib OPENDIA_WRITE_RAW=1
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/src/tests/class_tests/bin/OpenDIAScorer_arrays
P=/home/singjc/.jcode/scratch/e2e/peak
mkdir -p $P/score
for r in raw01 raw02; do
  for s in 0 1 2; do
    $B armE_$r $P/score/${r}_s$s.npy $s 8 armE_$r/heuristics.npy armE_$r/feat12.npy armE_$r/featD.npy cpp/featI_$r.npy \
      > $P/score/${r}_s$s.log 2>&1
    tail -1 $P/score/${r}_s$s.log
  done
done
