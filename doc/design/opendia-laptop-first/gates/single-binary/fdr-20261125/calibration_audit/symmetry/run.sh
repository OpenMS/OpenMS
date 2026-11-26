#!/usr/bin/env bash
set -euo pipefail
repo=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment
build=/home/singjc/Documents/github/OpenMS-wp03-release-build
scratch=/home/singjc/.jcode/scratch/symmetry-20261126
c++ -std=c++23 -O1 -fopenmp -I"$repo/src/topp/OpenDIA/include" -I"$repo/src/openms/include" -I"$build/src/openms/include" -I/usr/include/eigen3 "$scratch/symmetry.cpp" "$repo/src/topp/OpenDIA/source/OpenDIASearch.cpp" "$repo/src/topp/OpenDIA/source/OpenDIAPeakPicker.cpp" "$repo/src/topp/OpenDIA/source/OpenDIAScorer.cpp" -L"$build/lib" -Wl,-rpath,"$build/lib" -lOpenMS -o "$scratch/symmetry"
OMP_NUM_THREADS=2 "$scratch/symmetry" | tee "$scratch/results.txt"
