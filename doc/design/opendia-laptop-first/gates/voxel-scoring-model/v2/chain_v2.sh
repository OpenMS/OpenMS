#!/bin/bash
# PLAN_V2 chain: (export01 already running) -> export02 -> voxels 01/02 -> 4 arms per raw -> evaluate. STATUS_V2 = step.
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
V=$B/astral-voxel-pilot-20261110/v2
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
LIB=$B/astral-entrapment-20261107/lib/library.pqp
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
cd $V
echo wait-export01 > STATUS_V2
while [ "$(cat export01/STATUS 2>/dev/null)" != "done" ]; do sleep 20; done
grep -q IDENTICAL export01/CMP.txt || { echo export01-not-identical > STATUS_V2; exit 1; }
echo export02 > STATUS_V2
mkdir -p export02; cd export02
/usr/bin/time -v $BIN/OpenDIA -in /home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_02.raw \
  -lib $LIB -calibration $B/astral-entrapment-raw02-20261108/cal/cal.json -out search.tsv -features_out search.features.tsv \
  -min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10 \
  -out_xipm v2.xipm -xipm_precursors ../raw02/export_precursors.txt -xipm_mz_tolerance 20 -xipm_ms1_isotopes 3 > stdout.log 2> stderr.log
echo $? > EXIT
{ cmp search.features.tsv $B/astral-entrapment-raw02-20261108/search.features.tsv && echo IDENTICAL; } > CMP.txt 2>&1
cd $V
grep -q IDENTICAL export02/CMP.txt || { echo export02-not-identical > STATUS_V2; exit 1; }
for r in raw01 raw02; do
  echo voxels-$r > STATUS_V2
  /usr/bin/time -v $PY build_voxels.py $r > $r/voxels.log 2> $r/voxels.err || { echo voxels-$r-failed > STATUS_V2; exit 1; }
done
for r in raw01 raw02; do
  echo train-$r > STATUS_V2
  /usr/bin/time -v $PY train_v2.py $r V0 F0 S0 FS > $r/train.log 2> $r/train.err || { echo train-$r-failed > STATUS_V2; exit 1; }
done
echo evaluate > STATUS_V2
$PY evaluate_v2.py > EVALUATION_V2.txt 2> evaluate.err || { echo evaluate-failed > STATUS_V2; exit 1; }
echo done > STATUS_V2
