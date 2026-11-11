#!/bin/bash
# PLAN_V4 phase 2 data step: .xipm export of the 12 fragments for the v2 population precursors of each raw.
# The -lib is the 12-fragment population PQP; its search outputs are discarded (the population windows and their
# apexes stay those of the 6-fragment search). Same calibration and export options as v2/export0x.
# Waits for the v4 training chain (CPU and memory).
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
V=$B/astral-voxel-pilot-20261110/v2
L=$B/astral-lib12-20261111
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
cd $L; echo waiting > EXPORT_STATUS
while [ "$(cat $V/v4/STATUS 2>/dev/null)" != "done" ] && [ "$(cat $V/v4/STATUS 2>/dev/null | cut -c1-6)" != "failed" ]; do sleep 30; done
for r in 01 02; do
  CAL=$B/astral-design-raw01-20261107/window300k/cal.json; [ $r = 02 ] && CAL=$B/astral-entrapment-raw02-20261108/cal/cal.json
  mkdir -p export$r; cd export$r; echo export$r > ../EXPORT_STATUS
  /usr/bin/time -v $BIN/OpenDIA -in /home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_$r.raw \
    -lib $L/library12_pop.pqp -calibration $CAL -out search.tsv -features_out search.features.tsv \
    -min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10 \
    -out_xipm v12.xipm -xipm_precursors $V/raw$r/export_precursors.txt -xipm_mz_tolerance 20 -xipm_ms1_isotopes 3 \
    > stdout.log 2> stderr.log
  echo $? > EXIT; rm -f search.tsv search.features.tsv; cd $L
  [ "$(cat export$r/EXIT)" = 0 ] || { echo failed-export$r > EXPORT_STATUS; exit 1; }
done
echo done > EXPORT_STATUS
