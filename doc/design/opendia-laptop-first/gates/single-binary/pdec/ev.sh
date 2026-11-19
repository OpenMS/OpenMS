#!/bin/bash
# PRECHECK_EVIDENCE (diagnostic, seed 0 only): r~1.56 HeLa D with OpenDIA-predicted decoys (cached library),
# default vs -extra_scores. Entrapment evaluation afterwards with eval_pdec.py on the predicted library.
I=/home/singjc/.jcode/scratch/dia2/interf
W=/home/singjc/.jcode/scratch/dia2/gap7/ev
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIAWorkflow
D=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
L=/home/singjc/.jcode/scratch/dia2/ctest_cache/c/hela_entrap_codon.predicted_decoys.pqp
mkdir -p $W
/usr/bin/time -f "wall %E rss %M KB" $B -in $D -lib $L -out_dir $W/D_cur_s0 -threads 8 -seed 0 > $W/D_cur_s0.log 2>&1
echo "cur exit $? $(tail -1 $W/D_cur_s0.log)"
/usr/bin/time -f "wall %E rss %M KB" $B -in $D -lib $L -extra_scores -features_out -out_dir $W/D_new_s0 -threads 8 -seed 0 > $W/D_new_s0.log 2>&1
echo "new exit $? $(tail -1 $W/D_new_s0.log)"
grep -a "Scoring:" $W/D_cur_s0.log $W/D_new_s0.log | sed 's/rss_mb.*//'
