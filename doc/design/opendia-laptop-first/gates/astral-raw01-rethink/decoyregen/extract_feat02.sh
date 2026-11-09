#!/bin/bash
# Diagnostic re-extraction of raw 02 that KEEPS OSW's feature table (raw02_feat.osw), to localise which information our
# 20 subscores lack (gap_osw02.txt: 75% of the OSW human IDs we miss have the OSW peak as our rank-1 candidate but score
# a median 3.5 decoy-MAD below our FDP_c 1% cut). Same settings as extract.sh. OSW subscores are used only to diagnose
# which information is missing, never in a model we ship.
cd "$(dirname "$0")"; B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
R=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
mkdir -p tmp
/usr/bin/time -v systemd-run --user --scope -q -p MemoryMax=26G -p MemorySwapMax=0 \
$B/OpenSwathWorkflow -in $R/LFQ_Astral_DIA_15min_50ng_Human_02.raw -tr merged_cal02.pqp \
  -out_features raw02_feat.osw -readOptions onDisc -tempDirectory tmp -threads 12 \
  -enable_ms1 true -Scoring:stop_report_after_feature 5 \
  -Calibration:auto_irt:enabled false -Calibration:windows:estimate_rt false -Calibration:windows:estimate_mz false \
  -rt_extraction_window 200 -mz_extraction_window 34 -mz_extraction_window_ms1 55 \
  > extract_feat02.log 2>&1
echo $? > extract_feat02.EXIT
