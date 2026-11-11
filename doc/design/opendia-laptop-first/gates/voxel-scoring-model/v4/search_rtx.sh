#!/bin/bash
# A+B fix step 1: OpenDIA search with the extended calibration (extend_cal.py), same options as the base pipeline search
# (astral-design-raw01-20261107/pipeline/run.sh). No -out_xipm. Outputs in rtx/raw0X/.
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
V=$B/astral-voxel-pilot-20261110/v2
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
for r in 01 02; do
  CAL=$B/astral-design-raw01-20261107/window300k/cal.json; [ $r = 02 ] && CAL=$B/astral-entrapment-raw02-20261108/cal/cal.json
  O=$V/rtx/raw$r; mkdir -p $O; cd $O; echo search$r > $V/rtx/STATUS
  $PY $V/v4/extend_cal.py $CAL cal_ext.json > extend.log 2>&1 || { echo failed-cal > $V/rtx/STATUS; exit 1; }
  /usr/bin/time -v $BIN/OpenDIA -in /home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_$r.raw \
    -lib $B/astral-entrapment-20261107/lib/library.pqp -calibration cal_ext.json -out search.tsv -features_out search.features.tsv \
    -min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10 > search.stdout.log 2> search.stderr.log
  echo $? > EXIT; [ "$(cat EXIT)" = 0 ] || { echo failed-search$r > $V/rtx/STATUS; exit 1; }
done
echo done > $V/rtx/STATUS
