#!/bin/bash
# Raw 01 design diagnostic (see PLAN.md): -stage2_feature_set search+boundary. The CNN part is unchanged, so the logit,
# fold and window columns must be byte-identical to ../stage2-search-design/search (default-equivalence of the CNN path).
set -e
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
D=$B/astral-design-raw01-20261107
P=$D/pipeline
O=$D/stage2-byion-design
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
LIB=$B/astral-entrapment-20261107/lib/library.pqp
CAL=$D/window300k
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
cd $O
sha256sum $BIN/OpenDIACNN > binary.sha256
N=$(($(wc -l < $P/all_windows.tsv) - 1))
COMMON="-mode crossfit -in $IN -lib $LIB -windows $P/all_windows.tsv -seeds $P/seeds.tsv -max_samples $N -max_tensor_gib 20 -threads 5 -folds 3 -seed 42 -epochs 40 -patience 6 -batch_size 32 -paired_tdc -pooled_tdc -rt_calibration $CAL/cal.json -apply_calibration_coordinates -ss_iterations 0 -stage2"
mkdir -p byion
echo byion > STATUS
/usr/bin/time -v env PATH=/usr/bin:/bin $BIN/OpenDIACNN $COMMON -stage2_feature_set search+boundary+byion -stage2_search_features $P/search.features.tsv -model $O/byion/model -out $O/byion/crossfit.tsv -stage2_features_out $O/byion/features.tsv > byion/stdout.log 2> byion/stderr.log || { echo FAILED > STATUS; exit 1; }
{ cmp <(cut -f1-7,9 byion/crossfit.tsv) <(cut -f1-7,9 $D/stage2-boundary-design/boundary/crossfit.tsv) && echo "byion: sample/logit/fold/best/window columns IDENTICAL to stage2-boundary-design/boundary"; } > CMP.txt 2>&1 || true
{ cmp <(cut -f1-15 byion/features.tsv) <(cut -f1-15 $D/stage2-boundary-design/boundary/features.tsv) && echo "byion: first 15 feature columns (ids, cnn_z, search, boundary) IDENTICAL"; } >> CMP.txt 2>&1 || true
echo done > STATUS
