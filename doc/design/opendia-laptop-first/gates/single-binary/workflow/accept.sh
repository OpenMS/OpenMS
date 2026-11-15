#!/bin/bash
# PLAN_SINGLE_BINARY step 2-5 acceptance: OpenDIAWorkflow on the .raw files (arm-E calibration reused), seeds 0-4, both raws.
export DOTNET_ROOT=/home/singjc/.dotnet/runtime8 LD_LIBRARY_PATH=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib
G=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114
cd /home/singjc/.jcode/scratch/e2e
for r in 02 01; do
  for s in 0 1 2 3 4; do
    [ $r$s = 010 ] && continue
    R=$G/armE_raw$r
    echo "JCODE_PROGRESS {\"message\":\"raw$r seed $s\"}"
    /usr/bin/time -v bin/OpenDIAWorkflow -in /home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_$r.raw \
      -lib $R/library.pqp -scoring_lib $R/library12_pop.pqp -out_dir acc/raw${r}_s$s -threads 8 -seed $s \
      -calibration $R/cal/cal.json -anchors $R/cal/anchors.tsv > acc/raw${r}_s$s.log 2>&1 || echo "FAIL raw$r $s"
  done
done
echo DONE
