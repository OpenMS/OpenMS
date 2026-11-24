#!/usr/bin/env bash
set -euo pipefail
base=/home/singjc/.jcode/scratch/parquet/astral
work=/home/singjc/.jcode/scratch/compact-20261124
export DOTNET_ROOT=/home/singjc/.dotnet/runtime8
cp /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIA "$work/OpenDIA.release"
sha256sum "$work/OpenDIA.release" > "$work/release_binary.sha256"
raw=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
out="$work/out_release_cached_pqp"
mkdir "$out"
echo 'JCODE_PROGRESS {"message":"release binary raw01 seed0 automatic-cache PQP started"}'
/usr/bin/time -v -o "$work/release_cached_pqp.time.txt" "$work/OpenDIA.release" -in "$raw" -lib /home/singjc/.jcode/scratch/perf/lib/search.pqp -scoring_lib /home/singjc/.jcode/scratch/perf/lib/scoring.pqp -out_dir "$out" -threads 16 -seed 0 > "$work/release_cached_pqp.log" 2>&1
diff -rq "$base/out_base" "$out" | tee "$work/identity_release_cached_pqp.txt"
echo 'JCODE_PROGRESS {"message":"release automatic-cache outputs byte-identical to baseline"}'
