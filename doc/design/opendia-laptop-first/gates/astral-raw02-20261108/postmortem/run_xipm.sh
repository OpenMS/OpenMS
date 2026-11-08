#!/bin/bash
# Post-mortem only: raw 02 is spent. Re-runs the frozen OpenDIA search with -out_xipm on 290 selected precursors.
set -e
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_02.raw
LIB=../astral-entrapment-20261107/lib/library.pqp
echo running > STATUS
/usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration ../astral-entrapment-raw02-20261108/cal/cal.json -out search.tsv -min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10 -out_xipm selection.xipm -xipm_precursors xipm_precursors.txt > xipm.stdout.log 2> xipm.stderr.log
cmp search.tsv ../astral-entrapment-raw02-20261108/search.tsv && echo SEARCH_IDENTICAL > CMP || echo SEARCH_DIFFERS > CMP
echo done > STATUS
