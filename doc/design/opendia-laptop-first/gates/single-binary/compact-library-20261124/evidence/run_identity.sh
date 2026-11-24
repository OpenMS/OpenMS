#!/usr/bin/env bash
set -euo pipefail
base=/home/singjc/.jcode/scratch/parquet/astral
work=/home/singjc/.jcode/scratch/compact-20261124
mkdir -p "$work"
export DOTNET_ROOT=/home/singjc/.dotnet/runtime8
cp /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIA "$work/OpenDIA.compact"
sha256sum "$work/OpenDIA.compact" > "$work/binary.sha256"
raw=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
out="$work/out_parquet_attempt2"
mkdir "$out"
echo 'JCODE_PROGRESS {"message":"compact raw01 seed0 direct Parquet started"}'
/usr/bin/time -v -o "$work/parquet_attempt2.time.txt" "$work/OpenDIA.compact" -in "$raw" -lib /home/singjc/.jcode/scratch/perf/lib/search.opendia.parquet -scoring_lib /home/singjc/.jcode/scratch/perf/lib/scoring.opendia.parquet -out_dir "$out" -threads 16 -seed 0 > "$work/parquet_attempt2.log" 2>&1
diff -rq "$base/out_base" "$out" | tee "$work/identity_attempt2.txt"
echo 'JCODE_PROGRESS {"message":"compact outputs byte-identical to baseline"}'
