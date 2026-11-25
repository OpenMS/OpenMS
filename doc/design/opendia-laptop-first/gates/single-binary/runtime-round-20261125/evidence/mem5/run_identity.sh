#!/usr/bin/env bash
# Memory round 5: release extraction state before scoring, move features into the scorer.
# Same raw, seed, threads and libraries as the compact-20261124 identity arms; compare to the 3b1d90ed67 baseline.
set -euo pipefail
work=/home/singjc/.jcode/scratch/mem5-20261125
base=/home/singjc/.jcode/scratch/parquet/astral/out_base
mkdir -p "$work"
cp /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIA "$work/OpenDIA.mem5"
sha256sum "$work/OpenDIA.mem5" > "$work/binary.sha256"
git -C /home/singjc/Documents/github/OpenMS-r1-sparse-experiment diff > "$work/source.diff"
export DOTNET_ROOT=/home/singjc/.dotnet/runtime8
raw=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
out="$work/out_mem5"
mkdir "$out"
echo 'JCODE_PROGRESS {"message":"mem5 raw01 seed0 direct Parquet started"}'
/usr/bin/time -v -o "$work/mem5.time.txt" "$work/OpenDIA.mem5" -in "$raw" -lib /home/singjc/.jcode/scratch/perf/lib/search.opendia.parquet -scoring_lib /home/singjc/.jcode/scratch/perf/lib/scoring.opendia.parquet -out_dir "$out" -threads 16 -seed 0 > "$work/mem5.log" 2>&1
diff -rq "$base" "$out" | tee "$work/identity.txt"
test ! -s "$work/identity.txt"
echo 'JCODE_PROGRESS {"message":"mem5 outputs byte-identical to baseline"}'
