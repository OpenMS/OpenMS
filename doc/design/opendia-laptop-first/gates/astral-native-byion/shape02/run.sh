#!/bin/bash
set -e
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
R=$B/astral-entrapment-raw02-20261108
O=$B/astral-raw02-shape-diagnostic-20261109
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_02.raw
LIB=$B/astral-entrapment-20261107/lib/library.pqp
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
cd $O
sha256sum $BIN/OpenDIACNN > binary.sha256
N=$(($(wc -l < $R/all_windows.tsv) - 1))
COMMON="-mode crossfit -in $IN -lib $LIB -windows $R/all_windows.tsv -seeds $R/seeds.tsv -max_samples $N -max_tensor_gib 24 -threads 5 -folds 3 -seed 42 -epochs 40 -patience 6 -batch_size 32 -paired_tdc -pooled_tdc -rt_calibration $R/cal/cal.json -apply_calibration_coordinates -ss_iterations 0 -stage2 -stage2_search_features $R/search.features.tsv"
for arm in shape; do
  set_name=search+boundary; set_name=search+boundary+byion+shape
  mkdir -p $arm; echo $arm > STATUS
  /usr/bin/time -v env PATH=/usr/bin:/bin $BIN/OpenDIACNN $COMMON -stage2_feature_set $set_name -model $O/$arm/model -out $O/$arm/crossfit.tsv -stage2_features_out $O/$arm/features.tsv > $arm/stdout.log 2> $arm/stderr.log || { echo FAILED > STATUS; exit 1; }
done
echo done > STATUS
{ cmp <(cut -f1-16 shape/features.tsv) <(cut -f1-16 ../astral-raw02-byion-design-20261109/byion/features.tsv) && echo "feature cols 1-16 IDENTICAL to byion"; } > CMP.txt 2>&1 || true
