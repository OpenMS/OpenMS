#!/bin/bash
# Frozen Option B design, single evaluation on raw 02. Code: feature/opendia-clean-rewrite at 46be19510d.
# Runs once. Entrapment labels are joined only by evaluate.py after crossfit finishes.
set -e
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
O=$B/astral-entrapment-raw02-20261108
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
SRC=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment/src/topp/OpenDIA
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_02.raw
LIB=$B/astral-entrapment-20261107/lib/library.pqp
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
cd $O
[ -e STARTED ] && { echo "raw 02 evaluation already started once; refusing a second run" > REFUSED; exit 1; }
date -u +%FT%TZ > STARTED
mkdir -p cal crossfit
echo cal > STATUS
/usr/bin/time -v $BIN/OpenDIACalibrate -in $IN -lib $LIB -panel_targets 40000 -window_panel_targets 300000 -out cal/cal.json -out_anchors cal/anchors.tsv -out_diagnostics cal/diag.tsv > cal/stdout.log 2> cal/stderr.log
echo search > STATUS
/usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration cal/cal.json -out search.tsv -features_out search.features.tsv -min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10 > search.stdout.log 2> search.stderr.log
echo seeds > STATUS
python3 $SRC/prepare_trace_windows.py --library $LIB --anchors cal/anchors.tsv --candidates search.features.tsv --output seeds.tsv --all-windows all_windows.tsv --max-samples 6000 > seeds.log 2>&1
N=$(($(wc -l < all_windows.tsv) - 1))
echo cnn > STATUS
/usr/bin/time -v env PATH=/usr/bin:/bin $BIN/OpenDIACNN -mode crossfit -in $IN -lib $LIB -windows all_windows.tsv -seeds seeds.tsv -model $O/crossfit/model -out $O/crossfit/crossfit.tsv -max_samples $N -max_tensor_gib 24 -threads 5 -folds 3 -seed 42 -epochs 40 -patience 6 -batch_size 32 -paired_tdc -pooled_tdc -rt_calibration cal/cal.json -apply_calibration_coordinates -ss_iterations 0 > crossfit/stdout.log 2> crossfit/stderr.log
echo evaluate > STATUS
~/.venvs/redeem-external-dev-uv-20261022/bin/python $O/evaluate.py crossfit/crossfit.tsv $LIB > EVALUATION.txt
echo done > STATUS
