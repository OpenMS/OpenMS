#!/bin/bash
# PLAN_PAIRLIB data chain for one raw (run inside the SIF with the codon run root bound to the local path).
# Steps and options are those of xicids-20261112/chain.sh with the new option names of a8be5bed3c
# (-peak_picking xic removed, -xic_max_peaks -> -max_peaks) and the PLAN_PAIRLIB library.
# Usage: chain_raw.sh 01|02
set -o pipefail
r=$1
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
N=$B/pairlib-20261113
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
SRC=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment/src/topp/OpenDIA
LIB=$N/lib/library.pqp
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python
XIC="-min_transitions 1 -max_peaks 5"
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_$r.raw
st() { date -u +"%FT%TZ raw$r $1" >> $N/STATUS.log; }
fail() { st "FAIL $1"; exit 1; }
(cd $N && sha256sum -c FROZEN.sha256 > plan.check.raw$r 2>&1) || fail frozen-hash
R=$N/raw$r; mkdir -p $R/cal $R/pipeline $R/byion; cd $R
st cal-reused
# PLAN_DECOYGEN amendment 1: OpenDIACalibrate is not run. cal/cal.json and cal/anchors.tsv are the PAIRLIB run calibration
# (scripts/reuse_cal.py: cal.json unchanged, target anchors remapped by TRAML_ID). The calibration uses its own panel decoys and
# is a property of the run and the target library, which are unchanged; re-running it per arm failed or varied (RESULT notes).
[ -s cal/cal.json ] && [ -s cal/anchors.tsv ] || fail cal-missing
st search
/usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration cal/cal.json -out pipeline/search.tsv -features_out pipeline/search.features.tsv \
  $XIC > pipeline/search.stdout.log 2> pipeline/search.stderr.log || fail search
st seeds
$PY $SRC/prepare_trace_windows.py --library $LIB --anchors cal/anchors.tsv --candidates pipeline/search.features.tsv \
  --output pipeline/seeds.tsv --all-windows pipeline/all_windows.tsv --max-samples 6000 > pipeline/seeds.log 2>&1 || fail seeds
NW=$(($(wc -l < pipeline/all_windows.tsv) - 1))
st base
/usr/bin/time -v env PATH=/usr/bin:/bin $BIN/OpenDIACNN -mode crossfit -in $IN -lib $LIB -windows pipeline/all_windows.tsv -seeds pipeline/seeds.tsv \
  -max_samples $NW -max_tensor_gib 24 -threads 5 -folds 3 -seed 42 -epochs 40 -patience 6 -batch_size 32 -paired_tdc -pooled_tdc \
  -rt_calibration cal/cal.json -apply_calibration_coordinates -ss_iterations 0 -stage2 -stage2_feature_set search+boundary+byion \
  -stage2_search_features pipeline/search.features.tsv -model $R/byion/model -out $R/byion/crossfit.tsv -stage2_features_out $R/byion/features.tsv \
  > byion/stdout.log 2> byion/stderr.log || fail base
st population
$PY $N/scripts/make_population.py raw$r > population.log 2>&1 || fail population
# The 12-fragment PQP covers the populations of both raws; each raw job waits for the other's population, and raw01 builds it.
if [ $r = 01 ]; then
  while [ ! -s $N/raw02/population.tsv ] || [ ! -s $N/raw02/export_precursors.txt ]; do grep -q "raw02 FAIL" $N/STATUS.log && fail raw02-failed; sleep 30; done
  sleep 30; mkdir -p $N/lib12; st pqp12
  $PY $N/scripts/build_pqp12.py > $N/lib12/build_pqp12.log 2>&1 || fail pqp12
  touch $N/lib12/DONE
else
  while [ ! -e $N/lib12/DONE ]; do grep -q "raw01 FAIL" $N/STATUS.log && fail raw01-failed; sleep 30; done
fi
mkdir -p $N/export$r; cd $N/export$r
st export6
/usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration $R/cal/cal.json -out search.tsv -features_out search.features.tsv $XIC \
  -out_xipm v2.xipm -xipm_precursors $R/export_precursors.txt -xipm_mz_tolerance 20 -xipm_ms1_isotopes 3 > x6.stdout.log 2> x6.stderr.log || fail export6
cmp search.features.tsv $R/pipeline/search.features.tsv > CMP6.txt 2>&1 && echo IDENTICAL >> CMP6.txt || fail cmp6
rm -f search.tsv search.features.tsv
st export12
/usr/bin/time -v $BIN/OpenDIA -in $IN -lib $N/lib12/library12_pop.pqp -calibration $R/cal/cal.json -out search.tsv -features_out search.features.tsv $XIC \
  -out_xipm v12.xipm -xipm_precursors $R/export_precursors.txt -xipm_mz_tolerance 20 -xipm_ms1_isotopes 3 > x12.stdout.log 2> x12.stderr.log || fail export12
rm -f search.tsv search.features.tsv
cd $N
st voxels
$PY scripts/build_voxels.py raw$r > raw$r/voxels.log 2>&1 || fail voxels
st build12
$PY scripts/build12.py raw$r > raw$r/build12.log 2>&1 || fail build12
$PY scripts/features12.py raw$r > raw$r/features12.log 2>&1 || fail features12
st deficit
$PY scripts/deficit_features.py raw$r $N >> raw$r/deficit.log 2>&1 || fail deficit
rm -f export$r/v2.xipm export$r/v12.xipm
st data-done
