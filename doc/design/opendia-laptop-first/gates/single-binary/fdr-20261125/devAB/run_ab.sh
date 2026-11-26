#!/usr/bin/env bash
# FDR-20261125 arm A dev runs (Astral raw01/02 x seeds 0-4) plus a flag-off byte-identity check against raw01_s0 (OpenDIA.diag).
# One run at a time, 12 threads (19 GB free memory at launch; a run peaks near 13.6 GB). New paths only; nothing under fdr-20261125 is modified.
set -euo pipefail
H=/home/singjc/.jcode/scratch/fdr-ab; D=/home/singjc/.jcode/scratch/fdr-20261125; S=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment
mkdir -p $H/runs
[ -e $H/runs/OpenDIA.ab ] || cp /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIA $H/runs/OpenDIA.ab
sha256sum $H/runs/OpenDIA.ab > $H/runs/binary.sha256
git -C $S diff > $H/runs/source.diff
git -C $S rev-parse HEAD > $H/runs/source.head
export DOTNET_ROOT=/home/singjc/.dotnet/runtime8
one() {
  raw=$1; seed=$2; out=$3; shift 3
  [ -s "$out/precursors.tsv" ] && return 0
  mkdir -p "$out"
  /usr/bin/time -v -o "$out/time.txt" $H/runs/OpenDIA.ab -in "/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_${raw}.raw" \
    -lib /home/singjc/.jcode/scratch/perf/lib/search.opendia.parquet -scoring_lib /home/singjc/.jcode/scratch/perf/lib/scoring.opendia.parquet \
    -out_dir "$out" -threads 12 -seed "$seed" "$@" > "$out/run.log" 2>&1 || echo "FAIL $out"
  echo "JCODE_PROGRESS {\"message\":\"$(basename $out) done\"}"
}
one 01 0 $H/runs/off_raw01_s0
for f in precursors.tsv LFQ_Astral_DIA_15min_50ng_Human_01.precursors.tsv LFQ_Astral_DIA_15min_50ng_Human_01.calibration.json; do
  cmp "$D/raw01_s0/$f" "$H/runs/off_raw01_s0/$f" && echo "identical $f" || echo "DIFFERS $f"
done | tee $H/runs/identity.txt
for seed in 0 1 2 3 4; do for raw in 01 02; do
  one $raw $seed $H/runs/A_raw${raw}_s$seed -evidence_share -features_out
done; done
echo ALLDONE
