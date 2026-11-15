#!/bin/bash
# PLAN_CALREFINE2: CUR calibration (current code, raw02 anchors regenerated) and REF calibration reused, seeds 0-4, both raws.
export DOTNET_ROOT=/home/singjc/.dotnet/runtime8 LD_LIBRARY_PATH=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
RAW=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw
G=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114
E=/home/singjc/.jcode/scratch/e2e
cd $E && mkdir -p r2/cur02 r2/out
# Current-code calibration of raw02 with anchors (raw01 exists: cal_armE01).
if [ ! -s r2/cur02/anchors.tsv ]; then
  echo "JCODE_PROGRESS {\"message\":\"CUR raw02 calibration\"}"
  $BIN/OpenDIACalibrate -in $RAW/LFQ_Astral_DIA_15min_50ng_Human_02.raw -lib $G/armE_raw01/library.pqp -panel_targets 40000 \
    -window_panel_targets 300000 -out r2/cur02/cal.json -out_anchors r2/cur02/anchors.tsv > r2/cur02/log.txt 2>&1 || exit 1
fi
cmp r2/cur02/cal.json multi/out/LFQ_Astral_DIA_15min_50ng_Human_02.calibration.json && echo "CUR raw02 calibration reproduces the workflow run"
for s in 0 1 2 3 4; do
  for arm in CUR REF; do
    for r in 01 02; do
      if [ $arm = CUR ]; then
        if [ $r = 01 ]; then C=cal_armE01; else C=r2/cur02; fi
      else C=calrefine/armE_raw$r; fi
      O=r2/out/${arm}_raw${r}_s$s
      [ -s $O/LFQ_Astral_DIA_15min_50ng_Human_$r.precursors.tsv ] && continue
      echo "JCODE_PROGRESS {\"message\":\"$arm raw$r seed $s\"}"
      $BIN/OpenDIAWorkflow -in $RAW/LFQ_Astral_DIA_15min_50ng_Human_$r.raw -lib $G/armE_raw01/library.pqp -scoring_lib $G/armE_raw01/library12_pop.pqp \
        -out_dir $O -threads 8 -seed $s -calibration $C/cal.json -anchors $C/anchors.tsv > $O.log 2>&1 || echo "FAIL $arm $r $s"
    done
  done
done
echo DONE
