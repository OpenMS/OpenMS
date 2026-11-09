#!/bin/bash
# OSW extraction (XIC only, scores unused) of merged_cal{RUN}.pqp: targets + 3 decoy families, with library RT replaced
# by the run-calibrated RT that the reference analysis uses (calib_lib.py). OSW calibration is off, so the RT window
# (+-100 s) is centred there. 34 ppm MS2 / 55 ppm MS1, as in the reference OSW run. There is no m/z correction (the
# reference applied a quadratic correction before its 34 ppm window), which is a stated limit.
# An earlier attempt with OSW auto-iRT gave linear rsq 0.78 and off-centre windows (median |DIA-NN RT - centre| 50 s
# against 19 s for the reference XICs). Its result is kept as score02_offcentre.txt and is not used.
cd "$(dirname "$0")"; B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
R=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
mkdir -p tmp
for RUN in 02 01; do
  rm -f raw$RUN.osw raw$RUN.chrom.xic
  /usr/bin/time -v systemd-run --user --scope -q -p MemoryMax=26G -p MemorySwapMax=0 \
  $B/OpenSwathWorkflow -in $R/LFQ_Astral_DIA_15min_50ng_Human_$RUN.raw -tr merged_cal$RUN.pqp \
    -out_chrom raw$RUN.chrom.xic -out_features raw$RUN.osw -readOptions onDisc -tempDirectory tmp -threads 12 \
    -enable_ms1 true -Scoring:stop_report_after_feature 5 \
    -Calibration:auto_irt:enabled false -Calibration:windows:estimate_rt false -Calibration:windows:estimate_mz false \
    -rt_extraction_window 200 -mz_extraction_window 34 -mz_extraction_window_ms1 55 \
    > extract$RUN.log 2>&1
  echo $? > extract$RUN.EXIT
  rm -f raw$RUN.osw raw$RUN.osw-shm raw$RUN.osw-wal
done
