#!/usr/bin/env bash
# FDR diagnosis: default settings plus -features_out (adds network_score / network_q_value). Main outputs must stay
# byte-identical to the 3b1d90ed67 baseline on raw01 seed0.
set -euo pipefail
work=/home/singjc/.jcode/scratch/fdr-20261125
base=/home/singjc/.jcode/scratch/parquet/astral/out_base
cp /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIA "$work/OpenDIA.diag"
sha256sum "$work/OpenDIA.diag" > "$work/binary.sha256"
git -C /home/singjc/Documents/github/OpenMS-r1-sparse-experiment diff > "$work/source.diff"
export DOTNET_ROOT=/home/singjc/.dotnet/runtime8
for raw in 01 02; do
  out="$work/raw${raw}_s0"
  mkdir "$out"
  echo "JCODE_PROGRESS {\"message\":\"diag raw${raw} seed0 started\"}"
  /usr/bin/time -v -o "$out/time.txt" "$work/OpenDIA.diag" -in "/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_${raw}.raw" \
    -lib /home/singjc/.jcode/scratch/perf/lib/search.opendia.parquet -scoring_lib /home/singjc/.jcode/scratch/perf/lib/scoring.opendia.parquet \
    -out_dir "$out" -threads 16 -seed 0 -features_out > "$out/run.log" 2>&1
done
for f in precursors.tsv LFQ_Astral_DIA_15min_50ng_Human_01.precursors.tsv LFQ_Astral_DIA_15min_50ng_Human_01.calibration.json; do
  cmp "$base/$f" "$work/raw01_s0/$f"
done
echo 'JCODE_PROGRESS {"message":"diag raw01 main outputs byte-identical to baseline"}'
