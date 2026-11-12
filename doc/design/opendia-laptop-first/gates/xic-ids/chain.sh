#!/bin/bash
# PLAN_XIC_IDS data + score chain. Same steps and options as rtlin-20261111/chain_data.sh + chain_score.sh, except:
#  - binaries: commit 24dad4f842 (monotone RT map, fitted fragment tolerance, run-edge windows, peak_boundary block);
#  - every OpenDIA search (pipeline, 6- and 12-fragment exports) uses -peak_picking xic -xic_max_peaks 5 instead of
#    -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10 (boundary cap and min_profile_scans at defaults:
#    cap from cal.json, 2 scans).
# Scripts are copies of rtlin-20261111/scripts with the root repointed to xicids-20261112 (no other change).
set -o pipefail
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
N=$B/xicids-20261112
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
SRC=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment/src/topp/OpenDIA
LIB=$B/astral-entrapment-20261107/lib/library.pqp
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python
XIC="-min_transitions 1 -peak_picking xic -xic_max_peaks 5"
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
st() { echo "$1" > $N/STATUS; date -u +"%FT%TZ $1" >> $N/STATUS.log; }
fail() { st "failed-$1"; exit 1; }
cd $N
sha256sum -c PLAN_XIC_IDS.md.sha256 > PLAN_XIC_IDS.check 2>&1 || fail plan-changed
sha256sum $BIN/OpenDIA $BIN/OpenDIACalibrate $BIN/OpenDIACNN > binary.sha256
sha256sum -c binary.expected.sha256 > binary.check 2>&1 || fail binary-changed
for r in 01 02; do
  R=$N/raw$r; IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_$r.raw
  mkdir -p $R/cal $R/pipeline $R/byion; cd $R
  st cal$r
  /usr/bin/time -v $BIN/OpenDIACalibrate -in $IN -lib $LIB -panel_targets 40000 -window_panel_targets 300000 -out cal/cal.json \
    -out_anchors cal/anchors.tsv -out_diagnostics cal/diag.tsv > cal/stdout.log 2> cal/stderr.log || fail cal$r
  st search$r
  /usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration cal/cal.json -out pipeline/search.tsv -features_out pipeline/search.features.tsv \
    $XIC > pipeline/search.stdout.log 2> pipeline/search.stderr.log || fail search$r
  st seeds$r
  python3 $SRC/prepare_trace_windows.py --library $LIB --anchors cal/anchors.tsv --candidates pipeline/search.features.tsv \
    --output pipeline/seeds.tsv --all-windows pipeline/all_windows.tsv --max-samples 6000 > pipeline/seeds.log 2>&1 || fail seeds$r
  NW=$(($(wc -l < pipeline/all_windows.tsv) - 1))
  st base$r
  # Memory guard (no effect on results): BASE needed 8.9 KiB per window in the rtlin rebuild (23.1 GB / 2.60M windows).
  # Wait up to 3 h for MemAvailable >= 9.5 KiB x windows + 1 GB, then start; never kills other processes.
  NEED=$((NW * 95 / 10 + 1048576)); for i in $(seq 180); do AV=$(awk '/MemAvailable/{print $2}' /proc/meminfo); [ $AV -ge $NEED ] && break
    [ $i = 1 ] && date -u +"%FT%TZ base$r waiting for memory need ${NEED} kB available ${AV} kB" >> $N/STATUS.log; sleep 60; done
  [ $AV -ge $NEED ] || fail memory-base$r
  /usr/bin/time -v env PATH=/usr/bin:/bin $BIN/OpenDIACNN -mode crossfit -in $IN -lib $LIB -windows pipeline/all_windows.tsv -seeds pipeline/seeds.tsv \
    -max_samples $NW -max_tensor_gib 24 -threads 5 -folds 3 -seed 42 -epochs 40 -patience 6 -batch_size 32 -paired_tdc -pooled_tdc \
    -rt_calibration cal/cal.json -apply_calibration_coordinates -ss_iterations 0 -stage2 -stage2_feature_set search+boundary+byion \
    -stage2_search_features pipeline/search.features.tsv -model $R/byion/model -out $R/byion/crossfit.tsv -stage2_features_out $R/byion/features.tsv \
    > byion/stdout.log 2> byion/stderr.log || fail base$r
  st population$r
  $PY $N/scripts/make_population.py raw$r > population.log 2>&1 || fail population$r
done
cd $N; mkdir -p lib12
st pqp12
$PY scripts/build_pqp12.py > lib12/build_pqp12.log 2>&1 || fail pqp12
for r in 01 02; do
  R=$N/raw$r; IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_$r.raw
  mkdir -p $N/export$r; cd $N/export$r
  st export6_$r
  /usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration $R/cal/cal.json -out search.tsv -features_out search.features.tsv $XIC \
    -out_xipm v2.xipm -xipm_precursors $R/export_precursors.txt -xipm_mz_tolerance 20 -xipm_ms1_isotopes 3 > x6.stdout.log 2> x6.stderr.log || fail export6_$r
  cmp search.features.tsv $R/pipeline/search.features.tsv > CMP6.txt 2>&1 && echo IDENTICAL >> CMP6.txt || fail cmp6_$r
  rm -f search.tsv search.features.tsv
  st export12_$r
  /usr/bin/time -v $BIN/OpenDIA -in $IN -lib $N/lib12/library12_pop.pqp -calibration $R/cal/cal.json -out search.tsv -features_out search.features.tsv $XIC \
    -out_xipm v12.xipm -xipm_precursors $R/export_precursors.txt -xipm_mz_tolerance 20 -xipm_ms1_isotopes 3 > x12.stdout.log 2> x12.stderr.log || fail export12_$r
  rm -f search.tsv search.features.tsv
  cd $N
  st voxels$r
  $PY scripts/build_voxels.py raw$r > raw$r/voxels.log 2>&1 || fail voxels$r
  st build12_$r
  $PY scripts/build12.py raw$r > raw$r/build12.log 2>&1 || fail build12_$r
  $PY scripts/features12.py raw$r > raw$r/features12.log 2>&1 || fail features12_$r
  rm -f export$r/v2.xipm export$r/v12.xipm
done
st data-done
# Scoring: D11F -> W rule -> tail transform (unchanged recipe), then evaluation.
mkdir -p $N/score; cd $N
for r in raw01 raw02; do
  st train-$r
  TRAIN_THREADS=10 OMP_NUM_THREADS=10 /usr/bin/time -v $PY scripts/train_v4.py $r D11F > score/train_$r.log 2>&1 || fail train-$r
  $PY scripts/make_w.py $r D11F lin_D11F D11FW >> score/make_w.log 2>&1 || fail w-$r
  $PY scripts/calib_tail.py $r D11FW D11FWT >> score/calib_tail.log 2>&1 || fail t-$r
done
st evaluate
$PY scripts/evaluate_v2.py D11FW D11FWT > score/EVALUATION.txt 2> score/eval.err || fail eval
$PY scripts/null_auc.py BASE D11FWT > score/NULL.txt 2>> score/eval.err || fail null
$PY scripts/misses.py D11FWT > score/funnel_D11FWT.txt 2>> score/eval.err || fail funnel
$PY scripts/report_xic.py old=$B/astral-voxel-pilot-20261110/v2:D11FWT rtlin=$B/rtlin-20261111:D11FWT xic=$N:D11FWT \
  > score/REPORT_IDS.txt 2>> score/eval.err || fail report
st done
