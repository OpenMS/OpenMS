#!/bin/bash
# PLAN_CALREFINE criteria 1-2: OpenDIACalibrate (refined) on arm E raw01/02 and PAIRLIB raw01/02.
export DOTNET_ROOT=/home/singjc/.dotnet/runtime8 LD_LIBRARY_PATH=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIACalibrate
R=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw
L=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
cd /home/singjc/.jcode/scratch/e2e && mkdir -p calrefine
for job in armE:01:$L/gap-20261114/armE_raw01/library.pqp armE:02:$L/gap-20261114/armE_raw01/library.pqp pairlib:01:$L/pairlib-20261113/lib/library.pqp pairlib:02:$L/pairlib-20261113/lib/library.pqp; do
  IFS=: read name r lib <<< "$job"
  O=calrefine/${name}_raw$r; mkdir -p $O
  echo "JCODE_PROGRESS {\"message\":\"$name raw$r\"}"
  /usr/bin/time -v $B -in $R/LFQ_Astral_DIA_15min_50ng_Human_$r.raw -lib $lib -panel_targets 40000 -window_panel_targets 300000 \
    -out $O/cal.json -out_anchors $O/anchors.tsv > $O/log.txt 2>&1 || echo "FAIL $name $r"
done
echo DONE
