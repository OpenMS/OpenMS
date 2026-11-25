#!/usr/bin/env bash
# Runtime round 1: vectorised scorer forward pass (same per-unit summation order, bit-identical by construction).
# Same raw, seed, threads and libraries as the mem5 identity arm; compare to the 3b1d90ed67 baseline.
set -euo pipefail
work=/home/singjc/.jcode/scratch/rt1-20261125
base=/home/singjc/.jcode/scratch/parquet/astral/out_base
mkdir -p "$work"
cp /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIA "$work/OpenDIA.rt1"
sha256sum "$work/OpenDIA.rt1" > "$work/binary.sha256"
git -C /home/singjc/Documents/github/OpenMS-r1-sparse-experiment diff > "$work/source.diff"
export DOTNET_ROOT=/home/singjc/.dotnet/runtime8
raw=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
out="$work/out_rt1"
mkdir "$out"
echo 'JCODE_PROGRESS {"message":"rt1 raw01 seed0 direct Parquet started"}'
/usr/bin/time -v -o "$work/rt1.time.txt" "$work/OpenDIA.rt1" -in "$raw" -lib /home/singjc/.jcode/scratch/perf/lib/search.opendia.parquet -scoring_lib /home/singjc/.jcode/scratch/perf/lib/scoring.opendia.parquet -out_dir "$out" -threads 16 -seed 0 > "$work/rt1.log" 2>&1
diff -rq "$base" "$out" | tee "$work/identity.txt"
test ! -s "$work/identity.txt"
echo 'JCODE_PROGRESS {"message":"rt1 outputs byte-identical to baseline"}'
