#!/bin/bash
# Raw 01 design run only. Frozen search + -out_xipm for the matched-band selection.
set -e
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
LIB=../../astral-entrapment-20261107/lib/library.pqp
echo running > STATUS
/usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration ../window300k/cal.json -out search.tsv -min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10 -out_xipm selection.xipm -xipm_precursors precursors.txt > xipm.stdout.log 2> xipm.stderr.log
cmp search.tsv ../pipeline/search.tsv && echo SEARCH_IDENTICAL > CMP || echo SEARCH_DIFFERS > CMP
rm -f search.tsv
echo done > STATUS
