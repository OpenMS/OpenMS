#!/usr/bin/env bash
# More diagnostic seeds (same OpenDIA.diag binary, default options plus -features_out). Two runs in parallel, 10 threads each.
set -euo pipefail
work=/home/singjc/.jcode/scratch/fdr-20261125
sha256sum -c "$work/binary.sha256"
export DOTNET_ROOT=/home/singjc/.dotnet/runtime8
one() {
  raw=$1; seed=$2; out="$work/raw${raw}_s${seed}"
  mkdir "$out"
  /usr/bin/time -v -o "$out/time.txt" "$work/OpenDIA.diag" -in "/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_${raw}.raw" \
    -lib /home/singjc/.jcode/scratch/perf/lib/search.opendia.parquet -scoring_lib /home/singjc/.jcode/scratch/perf/lib/scoring.opendia.parquet \
    -out_dir "$out" -threads 10 -seed "$seed" -features_out > "$out/run.log" 2>&1
  echo "JCODE_PROGRESS {\"message\":\"raw${raw} s${seed} done\"}"
}
for seed in 1 2 3 4; do
  one 01 $seed & one 02 $seed & wait
done
