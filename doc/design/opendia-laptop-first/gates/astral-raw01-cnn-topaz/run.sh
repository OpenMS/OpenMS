#!/bin/bash
# Raw 01 TOPAZ design runs (PLAN.md): A = -cnn_arch topaz (seeds), B = + -cnn_training mil.
# Base: the frozen stage-2 search command (../stage2-search-design/run.sh, search arm).
set -e
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
D=$B/astral-design-raw01-20261107
P=$D/pipeline
O=$D/cnn-topaz-design
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
LIB=$B/astral-entrapment-20261107/lib/library.pqp
CAL=$D/window300k
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
cd $O
N=$(($(wc -l < $P/all_windows.tsv) - 1))
COMMON="-mode crossfit -in $IN -lib $LIB -windows $P/all_windows.tsv -seeds $P/seeds.tsv -max_samples $N -max_tensor_gib 20 -threads 5 -folds 3 -seed 42 -epochs 40 -patience 6 -batch_size 32 -paired_tdc -pooled_tdc -rt_calibration $CAL/cal.json -apply_calibration_coordinates -ss_iterations 0 -stage2 -stage2_feature_set search -stage2_search_features $P/search.features.tsv"
for arm in ${ARMS:-A B}; do
  case $arm in
    A) EXTRA="-cnn_arch topaz" ;;
    B) EXTRA="-cnn_arch topaz -cnn_training mil" ;;
  esac
  mkdir -p $arm
  echo $arm > STATUS
  /usr/bin/time -v env PATH=/usr/bin:/bin $BIN/OpenDIACNN $COMMON $EXTRA -model $O/$arm/model -out $O/$arm/crossfit.tsv -stage2_features_out $O/$arm/features.tsv > $arm/stdout.log 2> $arm/stderr.log || { echo FAILED_$arm > STATUS; exit 1; }
done
echo done > STATUS
