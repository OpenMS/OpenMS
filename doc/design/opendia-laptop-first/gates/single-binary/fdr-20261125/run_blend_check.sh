#!/usr/bin/env bash
# FDR-20261125 blend: (1) default byte-identity vs the diag run raw01_s0, (2) -gbm_blend 0.5 on raw01_s0 (offline sim: N 2047, E 11 approx).
set -euo pipefail
work=/home/singjc/.jcode/scratch/fdr-20261125
cp /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIA "$work/OpenDIA.blend"
sha256sum "$work/OpenDIA.blend" > "$work/blend_binary.sha256"
git -C /home/singjc/Documents/github/OpenMS-r1-sparse-experiment diff > "$work/blend_source.diff"
export DOTNET_ROOT=/home/singjc/.dotnet/runtime8
run() {
  out=$1; shift
  mkdir "$out"
  "$work/OpenDIA.blend" -in /home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw \
    -lib /home/singjc/.jcode/scratch/perf/lib/search.opendia.parquet -scoring_lib /home/singjc/.jcode/scratch/perf/lib/scoring.opendia.parquet \
    -out_dir "$out" -threads 16 -seed 0 "$@" > "$out/run.log" 2>&1
}
run "$work/blend_off_raw01_s0"
echo 'JCODE_PROGRESS {"message":"flag-off run done"}'
for f in precursors.tsv LFQ_Astral_DIA_15min_50ng_Human_01.precursors.tsv LFQ_Astral_DIA_15min_50ng_Human_01.calibration.json; do
  cmp "$work/raw01_s0/$f" "$work/blend_off_raw01_s0/$f" && echo "identical $f"
done
run "$work/blend_w05_raw01_s0" -gbm_blend 0.5 -features_out
echo 'JCODE_PROGRESS {"message":"blend 0.5 run done"}'
grep -h "GBM blend\|q <= 0.01" "$work/blend_w05_raw01_s0/run.log"
