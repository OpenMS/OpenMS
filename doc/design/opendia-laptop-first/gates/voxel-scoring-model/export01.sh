#!/bin/bash
# Voxel pilot raw 01 export: raw MS1/MS2 points for the training and evaluation precursors (export_precursors.txt).
# Same search options as the base pipeline, so the -features_out candidate rows keep the base sample ids (checked by
# comparing with the base search.features.tsv). Waits for the round 7 chain to finish (memory).
set -e
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
D=$B/astral-design-raw01-20261107
O=$B/astral-voxel-pilot-20261110/export01
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
LIB=$B/astral-entrapment-20261107/lib/library.pqp
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
mkdir -p $O && cd $O
echo waiting > STATUS
while [ "$(cat $B/astral-single-scorer-20261109/STATUS7 2>/dev/null)" != "done" ]; do sleep 20; done
echo export > STATUS
/usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration $D/window300k/cal.json -out search.tsv -features_out search.features.tsv \
  -min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10 \
  -out_xipm pilot.xipm -xipm_precursors ../export_precursors.txt -xipm_mz_tolerance 20 -xipm_ms1_isotopes 3 \
  > stdout.log 2> stderr.log
echo $? > EXIT
cmp search.features.tsv $D/pipeline/search.features.tsv > CMP.txt 2>&1 && echo IDENTICAL >> CMP.txt || echo DIFFERENT >> CMP.txt
echo done > STATUS
