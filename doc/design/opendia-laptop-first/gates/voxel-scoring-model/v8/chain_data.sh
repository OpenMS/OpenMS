#!/bin/bash
# rtlin rebuild: the full raw 01/02 pipeline with the linear-end RTMap (OpenMS-r1-sparse-experiment working tree,
# binaries in OpenMS-r1-full-ninja-build/bin built 2026-11-11 15:3x). Same options as the frozen pipeline
# (astral-entrapment-raw02-20261108/run.sh, astral-design-raw01-20261107/{window300k,pipeline,stage2-byion-design}/run.sh,
# astral-voxel-pilot-20261110/v2/export01.sh, astral-lib12-20261111/export12.sh).
# Steps per raw: calibrate -> search -> seeds/all_windows -> BASE (OpenDIACNN -stage2 byion) -> population.
# Then both raws: 12-fragment PQP for the new populations -> per raw: 6- and 12-fragment xipm export, arrays, feat12.
# Scoring (D11F, W rule, tail transform) and evaluation run only after PLAN_V8 is frozen (chain_score.sh).
set -o pipefail
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
N=$B/rtlin-20261111
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
SRC=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment/src/topp/OpenDIA
LIB=$B/astral-entrapment-20261107/lib/library.pqp
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
st() { echo "$1" > $N/STATUS; date -u +"%FT%TZ $1" >> $N/STATUS.log; }
fail() { st "failed-$1"; exit 1; }
cd $N
sha256sum $BIN/OpenDIA $BIN/OpenDIACalibrate $BIN/OpenDIACNN > binary.sha256
for r in 01 02; do
  R=$N/raw$r; IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_$r.raw
  mkdir -p $R/cal $R/pipeline/crossfit $R/byion; cd $R
  st cal$r
  /usr/bin/time -v $BIN/OpenDIACalibrate -in $IN -lib $LIB -panel_targets 40000 -window_panel_targets 300000 -out cal/cal.json \
    -out_anchors cal/anchors.tsv -out_diagnostics cal/diag.tsv > cal/stdout.log 2> cal/stderr.log || fail cal$r
  st search$r
  /usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration cal/cal.json -out pipeline/search.tsv -features_out pipeline/search.features.tsv \
    -min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10 > pipeline/search.stdout.log 2> pipeline/search.stderr.log || fail search$r
  st seeds$r
  python3 $SRC/prepare_trace_windows.py --library $LIB --anchors cal/anchors.tsv --candidates pipeline/search.features.tsv \
    --output pipeline/seeds.tsv --all-windows pipeline/all_windows.tsv --max-samples 6000 > pipeline/seeds.log 2>&1 || fail seeds$r
  NW=$(($(wc -l < pipeline/all_windows.tsv) - 1))
  st base$r
  /usr/bin/time -v env PATH=/usr/bin:/bin $BIN/OpenDIACNN -mode crossfit -in $IN -lib $LIB -windows pipeline/all_windows.tsv -seeds pipeline/seeds.tsv \
    -max_samples $NW -max_tensor_gib 24 -threads 5 -folds 3 -seed 42 -epochs 40 -patience 6 -batch_size 32 -paired_tdc -pooled_tdc \
    -rt_calibration cal/cal.json -apply_calibration_coordinates -ss_iterations 0 -stage2 -stage2_feature_set search+boundary+byion \
    -stage2_search_features pipeline/search.features.tsv -model $R/byion/model -out $R/byion/crossfit.tsv -stage2_features_out $R/byion/features.tsv \
    > byion/stdout.log 2> byion/stderr.log || fail base$r
  st population$r
  $PY $N/scripts/make_population.py raw$r > population.log 2>&1 || fail population$r
done
cd $N; mkdir -p lib12 v3seq v4
st pqp12
$PY scripts/build_pqp12.py > lib12/build_pqp12.log 2>&1 || fail pqp12
for r in 01 02; do
  R=$N/raw$r; IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_$r.raw
  mkdir -p $N/export$r; cd $N/export$r
  st export6_$r
  /usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration $R/cal/cal.json -out search.tsv -features_out search.features.tsv \
    -min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10 \
    -out_xipm v2.xipm -xipm_precursors $R/export_precursors.txt -xipm_mz_tolerance 20 -xipm_ms1_isotopes 3 > x6.stdout.log 2> x6.stderr.log || fail export6_$r
  cmp search.features.tsv $R/pipeline/search.features.tsv > CMP6.txt 2>&1 && echo IDENTICAL >> CMP6.txt || fail cmp6_$r
  rm -f search.tsv search.features.tsv
  st export12_$r
  /usr/bin/time -v $BIN/OpenDIA -in $IN -lib $N/lib12/library12_pop.pqp -calibration $R/cal/cal.json -out search.tsv -features_out search.features.tsv \
    -min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10 \
    -out_xipm v12.xipm -xipm_precursors $R/export_precursors.txt -xipm_mz_tolerance 20 -xipm_ms1_isotopes 3 > x12.stdout.log 2> x12.stderr.log || fail export12_$r
  rm -f search.tsv search.features.tsv
  cd $N
  st voxels$r
  $PY scripts/build_voxels.py raw$r > raw$r/voxels.log 2>&1 || fail voxels$r
  st build12_$r
  $PY scripts/build12.py raw$r > raw$r/build12.log 2>&1 || fail build12_$r
  $PY scripts/features12.py raw$r > raw$r/features12.log 2>&1 || fail features12_$r
  rm -f export$r/v2.xipm export$r/v12.xipm  # arrays built; xipm exports are large (about 1.5 GB each)
done
st data-done
