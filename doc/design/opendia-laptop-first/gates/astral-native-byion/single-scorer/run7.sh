#!/bin/bash
# PLAN7.md chain: HeLa default equivalence (fresh dir), then raw01/P1, raw02/P1 (sequential).
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
O=$B/astral-single-scorer-20261109
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
LIB=$B/astral-entrapment-20261107/lib/library.pqp
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
cd $O
sha256sum $BIN/OpenDIACNN > binary7.sha256
H=$O/hela-default-cmp7
mkdir -p $H
awk -v h="$H" '/^O=/{print "O=" h; next} {print}' $B/astral-design-raw01-20261107/cnn-topaz-design/hela-default-cmp4/run.sh > $H/run.sh
echo hela > STATUS7
bash $H/run.sh
grep -q IDENTICAL $H/CMP.txt || { echo FAILED_hela_equivalence > STATUS7; exit 1; }
COMMON="-mode crossfit -lib $LIB -max_tensor_gib 24 -threads 5 -folds 3 -seed 42 -epochs 40 -patience 6 -batch_size 32 -paired_tdc -pooled_tdc -apply_calibration_coordinates -ss_iterations 0 -stage2 -stage2_feature_set search+boundary+byion -stage2_refined_picking_min_r 0.75"
for raw in raw01 raw02; do
  if [ $raw = raw01 ]; then
    P=$B/astral-design-raw01-20261107/pipeline; CAL=$B/astral-design-raw01-20261107/window300k/cal.json
    IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
    BASE=$B/astral-design-raw01-20261107/stage2-byion-design/byion
  else
    P=$B/astral-entrapment-raw02-20261108; CAL=$P/cal/cal.json
    IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_02.raw
    BASE=$B/astral-raw02-byion-design-20261109/byion
  fi
  N=$(($(wc -l < $P/all_windows.tsv) - 1))
  W=$O/$raw/P1; mkdir -p $W
  echo $raw/P1 > STATUS7
  /usr/bin/time -v env PATH=/usr/bin:/bin $BIN/OpenDIACNN $COMMON -in $IN -windows $P/all_windows.tsv -seeds $P/seeds.tsv -max_samples $N \
    -rt_calibration $CAL -stage2_search_features $P/search.features.tsv -model $W/model -out $W/crossfit.tsv -stage2_features_out $W/features.tsv \
    -quant_out $W/quant.tsv > $W/stdout.log 2> $W/stderr.log || { echo FAILED_$raw > STATUS7; exit 1; }
  { cmp <(cut -f1-6,9 $W/crossfit.tsv) <(cut -f1-6,9 $BASE/crossfit.tsv) && echo "$raw/P1: crossfit cols 1-6,9 IDENTICAL to base"; } >> CMP7.txt 2>&1
  { cmp <(cut -f1-9,15-16 $W/features.tsv) <(cut -f1-9,15-16 $BASE/features.tsv) && echo "$raw/P1: features cols 1-9,15-16 IDENTICAL to base"; } >> CMP7.txt 2>&1
done
echo done > STATUS7
