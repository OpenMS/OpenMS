#!/usr/bin/env bash
# Build the A/B replay harness against the current OpenDIA_impl (same compile flags as the tool's sources).
set -euo pipefail
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build
S=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment/src/topp/OpenDIA
H=/home/singjc/.jcode/scratch/fdr-ab
g++ -O3 -DNDEBUG -std=gnu++23 -mssse3 -ffp-contract=off -fopenmp -I$S/include -o $H/gbm_replay $H/gbm_replay.cpp \
  $B/lib/libOpenDIA_impl.a -Wl,-rpath,$B/lib -L/usr/lib/gcc/x86_64-linux-gnu/14 /usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so -lpthread
echo built
