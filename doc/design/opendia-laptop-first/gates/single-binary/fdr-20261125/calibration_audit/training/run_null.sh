#!/usr/bin/env bash
set -euo pipefail
repo=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment/src/topp/OpenDIA
here=$(cd "$(dirname "$0")" && pwd)
c++ -O2 -std=c++17 -fopenmp -I"$repo/include" "$here/null_stress.cpp" "$repo/source/OpenDIAGBM.cpp" "$repo/source/OpenDIAScorer.cpp" -o "$here/null_stress"
printf 'seed\tarm\taccepted_targets\tknown_null_targets\taccepted_decoys\tFDP\n' > "$here/null_stress.tsv"
: > "$here/null_stress_failures.txt"
for seed in {0..9}; do
  if "$here/null_stress" "$seed" > "$here/null_${seed}.tsv" 2> "$here/null_${seed}.err"; then
    cat "$here/null_${seed}.tsv" >> "$here/null_stress.tsv"
  else
    printf '%s FAILED\n' "$seed" >> "$here/null_stress_failures.txt"
  fi
done
cat "$here/null_stress.tsv" "$here/null_stress_failures.txt"
