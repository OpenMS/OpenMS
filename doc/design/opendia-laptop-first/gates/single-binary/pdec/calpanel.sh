#!/bin/bash
# PRECHECK_CALPANEL (diagnostic): does a larger calibration panel restore anchors on r~7 D?
# 40000 is the workflow's panel (must reproduce 140 used anchors). Local, 8 threads.
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIACalibrate
L=/home/singjc/.jcode/scratch/dia2/entlib10/cut/hela_r5.predicted_decoys.pqp
D=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
W=/home/singjc/.jcode/scratch/dia2/gap7/cal
mkdir -p $W
for P in 40000 160000 320000; do
  /usr/bin/time -f "wall %E rss %M KB" $B -in $D -lib $L -panel_targets $P -out $W/p$P.json -out_anchors $W/p$P.anchors.tsv -threads 8 > $W/p$P.log 2>&1
  echo "panel $P exit $? $(tail -1 $W/p$P.log)"
  grep -o '"stage1": {[^}]*}' $W/p$P.json | cut -c1-120
  grep -o '"stage3": {[^}]*}' $W/p$P.json
  grep -o '"anchors_used": [0-9]*\|"tolerance_half_width_s": [0-9.]*\|"abs_residual_p95_s": [0-9.]*' $W/p$P.json | tr '\n' ' '; echo
done
