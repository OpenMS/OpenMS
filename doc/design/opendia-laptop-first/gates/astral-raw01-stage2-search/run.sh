#!/bin/bash
# 1) -stage2 default (extraction set) with the refactored binary: must be byte-identical to ../stage2-design.
# 2) -stage2_feature_set search (PLAN.md).
set -e
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
D=$B/astral-design-raw01-20261107
P=$D/pipeline
O=$D/stage2-search-design
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
LIB=$B/astral-entrapment-20261107/lib/library.pqp
CAL=$D/window300k
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
cd $O
N=$(($(wc -l < $P/all_windows.tsv) - 1))
COMMON="-mode crossfit -in $IN -lib $LIB -windows $P/all_windows.tsv -seeds $P/seeds.tsv -max_samples $N -max_tensor_gib 20 -threads 5 -folds 3 -seed 42 -epochs 40 -patience 6 -batch_size 32 -paired_tdc -pooled_tdc -rt_calibration $CAL/cal.json -apply_calibration_coordinates -ss_iterations 0 -stage2"
mkdir -p default search
echo default > STATUS
/usr/bin/time -v env PATH=/usr/bin:/bin $BIN/OpenDIACNN $COMMON -model $O/default/model -out $O/default/crossfit.tsv -stage2_features_out $O/default/features.tsv > default/stdout.log 2> default/stderr.log || { echo FAILED_default > STATUS; exit 1; }
{ cmp default/crossfit.tsv $D/stage2-design/crossfit.tsv && cmp default/features.tsv $D/stage2-design/features.tsv && echo "default -stage2 output IDENTICAL to stage2-design"; } > CMP_default.txt 2>&1 || true
echo search > STATUS
/usr/bin/time -v env PATH=/usr/bin:/bin $BIN/OpenDIACNN $COMMON -stage2_feature_set search -stage2_search_features $P/search.features.tsv -model $O/search/model -out $O/search/crossfit.tsv -stage2_features_out $O/search/features.tsv > search/stdout.log 2> search/stderr.log || { echo FAILED_search > STATUS; exit 1; }
{ cmp <(cut -f1-6,9 search/crossfit.tsv) <(cut -f1-6,9 $D/stage2-design/crossfit.tsv) && echo "search: logit/fold/window columns IDENTICAL"; } > CMP_search.txt 2>&1 || true
echo done > STATUS
