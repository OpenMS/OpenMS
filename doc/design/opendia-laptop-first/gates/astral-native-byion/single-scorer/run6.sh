#!/bin/bash
# PLAN6.md chain: wait for the HeLa default equivalence (hela-default-cmp6), then raw01/R1, raw02/R1, raw01/R2, raw02/R2 (sequential).
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
O=$B/astral-single-scorer-20261109
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
LIB=$B/astral-entrapment-20261107/lib/library.pqp
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
cd $O
sha256sum $BIN/OpenDIACNN > binary6.sha256
H=$O/hela-default-cmp6
echo hela > STATUS6
while [ ! -f $H/DONE ]; do sleep 20; done
grep -q IDENTICAL $H/CMP.txt || { echo FAILED_hela_equivalence > STATUS6; exit 1; }
cmp -s $H/binary.sha256 binary6.sha256 || { echo FAILED_binary_changed > STATUS6; exit 1; }
COMMON="-mode crossfit -lib $LIB -max_tensor_gib 24 -threads 5 -folds 3 -seed 42 -epochs 40 -patience 6 -batch_size 32 -paired_tdc -pooled_tdc -apply_calibration_coordinates -ss_iterations 0 -stage2 -stage2_feature_set search+boundary+byion+refine"
for job in raw01/R1 raw02/R1 raw01/R2 raw02/R2; do
  raw=${job%/*}; arm=${job#*/}
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
  EXTRA=""; [ $arm = R2 ] && EXTRA="-trace_refinement_min_r 0.75"
  W=$O/$raw/$arm; mkdir -p $W
  echo $job > STATUS6
  /usr/bin/time -v env PATH=/usr/bin:/bin $BIN/OpenDIACNN $COMMON $EXTRA -in $IN -windows $P/all_windows.tsv -seeds $P/seeds.tsv -max_samples $N \
    -rt_calibration $CAL -stage2_search_features $P/search.features.tsv -model $W/model -out $W/crossfit.tsv -stage2_features_out $W/features.tsv \
    > $W/stdout.log 2> $W/stderr.log || { echo FAILED_$job > STATUS6; exit 1; }
  if [ $arm = R1 ]; then
    { cmp <(cut -f1-7,9 $W/crossfit.tsv) <(cut -f1-7,9 $BASE/crossfit.tsv) && echo "$job: crossfit cols 1-7,9 IDENTICAL to base"; } >> CMP6.txt 2>&1
    { cmp <(cut -f1-16 $W/features.tsv) <(cut -f1-16 $BASE/features.tsv) && echo "$job: features cols 1-16 IDENTICAL to base"; } >> CMP6.txt 2>&1
  fi
done
echo done > STATUS6
