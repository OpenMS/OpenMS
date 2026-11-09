#!/bin/bash
# PLAN.md chain: HeLa default equivalence, then raw01 S0, S1, raw02 S0, S1 (sequential, about 25 GB each).
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
O=$B/astral-single-scorer-20261109
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
LIB=$B/astral-entrapment-20261107/lib/library.pqp
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
cd $O
sha256sum $BIN/OpenDIACNN > binary.sha256
# Default equivalence (frozen HeLa command, byte-identical crossfit.tsv required).
H=$O/hela-default-cmp
mkdir -p $H
sed "s#^O=.*#O=$H#" $B/astral-design-raw01-20261107/cnn-topaz-design/hela-default-cmp4/run.sh > $H/run.sh
echo hela > STATUS
bash $H/run.sh
grep -q IDENTICAL $H/CMP.txt || { echo FAILED_hela_equivalence > STATUS; exit 1; }
COMMON="-mode crossfit -lib $LIB -max_tensor_gib 24 -threads 5 -folds 3 -seed 42 -epochs 40 -patience 6 -batch_size 32 -paired_tdc -pooled_tdc -apply_calibration_coordinates -cnn_arch topaz -topaz_width 8 -weight_decay 0.01 -topaz_features search+boundary+byion -ss_iterations 2 -ss_fdr 0.01"
for raw in raw01 raw02; do
  if [ $raw = raw01 ]; then
    P=$B/astral-design-raw01-20261107/pipeline; CAL=$B/astral-design-raw01-20261107/window300k/cal.json
    IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
  else
    P=$B/astral-entrapment-raw02-20261108; CAL=$P/cal/cal.json
    IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_02.raw
  fi
  N=$(($(wc -l < $P/all_windows.tsv) - 1))
  for arm in S0 S1; do
    EXTRA=""; [ $arm = S1 ] && EXTRA="-augment_probability 0.5 -augment_scale_min 0.1"
    W=$O/$raw/$arm; mkdir -p $W
    echo $raw/$arm > STATUS
    /usr/bin/time -v env PATH=/usr/bin:/bin $BIN/OpenDIACNN $COMMON $EXTRA -in $IN -windows $P/all_windows.tsv -seeds $P/seeds.tsv -max_samples $N \
      -rt_calibration $CAL -stage2_search_features $P/search.features.tsv -model $W/model -out $W/crossfit.tsv > $W/stdout.log 2> $W/stderr.log \
      || { echo FAILED_$raw/$arm > STATUS; exit 1; }
  done
done
echo done > STATUS
