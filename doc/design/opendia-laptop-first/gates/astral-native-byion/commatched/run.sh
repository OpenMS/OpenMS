#!/bin/bash
# k11 native runs (PLAN.md): only the decoys differ. Calibration (cal.json) is reused; anchors are remapped by TRAML_ID.
set -e
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
O=$B/astral-commatched-decoys-20261109
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
SRC=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment/src/topp/OpenDIA
LIB=$O/library.pqp
PY=$1
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
cd $O
sha256sum $BIN/OpenDIACNN $BIN/OpenDIA > binary.sha256
for run in raw01 raw02; do
  if [ $run = raw01 ]; then
    IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
    CAL=$B/astral-design-raw01-20261107/window300k; TG=20
  else
    IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_02.raw
    CAL=$B/astral-entrapment-raw02-20261108/cal; TG=24
  fi
  mkdir -p $run/crossfit; echo $run search > STATUS
  $PY $O/remap_anchors.py $B/astral-entrapment-20261107/lib/library.pqp $LIB $CAL/anchors.tsv $O/$run/anchors.tsv > $O/$run/remap.log
  /usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration $CAL/cal.json -out $O/$run/search.tsv -features_out $O/$run/search.features.tsv -min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10 > $O/$run/search.stdout.log 2> $O/$run/search.stderr.log
  echo $run seeds > STATUS
  python3 $SRC/prepare_trace_windows.py --library $LIB --anchors $O/$run/anchors.tsv --candidates $O/$run/search.features.tsv --output $O/$run/seeds.tsv --all-windows $O/$run/all_windows.tsv --max-samples 6000 > $O/$run/seeds.log 2>&1
  N=$(($(wc -l < $O/$run/all_windows.tsv) - 1))
  echo $run cnn > STATUS
  /usr/bin/time -v env PATH=/usr/bin:/bin $BIN/OpenDIACNN -mode crossfit -in $IN -lib $LIB -windows $O/$run/all_windows.tsv -seeds $O/$run/seeds.tsv -max_samples $N -max_tensor_gib $TG -threads 5 -folds 3 -seed 42 -epochs 40 -patience 6 -batch_size 32 -paired_tdc -pooled_tdc -rt_calibration $CAL/cal.json -apply_calibration_coordinates -ss_iterations 0 -stage2 -stage2_feature_set search+boundary+byion -stage2_search_features $O/$run/search.features.tsv -model $O/$run/crossfit/model -out $O/$run/crossfit/crossfit.tsv -stage2_features_out $O/$run/crossfit/features.tsv > $O/$run/crossfit/stdout.log 2> $O/$run/crossfit/stderr.log
done
echo done > STATUS
