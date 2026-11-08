#!/bin/bash
# Fragment ablation: identical search as pipeline/run.sh, library reduced to the top-4 transitions per precursor.
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107; D=$B/astral-design-raw01-20261107
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet; export PATH="$DOTNET_ROOT:$PATH"
cd $D/rethink/frag4
/usr/bin/time -v $BIN/OpenDIA -in $IN -lib lib4.pqp -calibration $D/window300k/cal.json -out search.tsv -features_out search.features.tsv -min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10 > search.stdout.log 2> search.stderr.log
echo $? > EXIT
