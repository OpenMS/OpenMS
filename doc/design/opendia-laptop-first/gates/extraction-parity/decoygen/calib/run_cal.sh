#!/bin/bash
# Calibration probe: OpenDIACalibrate on one library and raw (same options as chain_raw.sh). Usage: run_cal.sh <outdir> <library.pqp> <raw 01|02>
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_$3.raw
mkdir -p "$1" && cd "$1" || exit 1
/usr/bin/time -v "$B/OpenDIACalibrate" -in "$IN" -lib "$2" -panel_targets 40000 -window_panel_targets 300000 -out cal.json \
  -out_anchors anchors.tsv -out_diagnostics diag.tsv > stdout.log 2> stderr.log
echo $? > EXIT
