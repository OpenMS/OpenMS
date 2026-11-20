#!/bin/bash
# PRECHECK_IMNARROW (diagnostic, seed 0): r~1.56 HeLa D, default vs narrower calibrated 1/K0 window (search, grid, interference).
W=/home/singjc/.jcode/scratch/dia2/gap7/ev
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIAWorkflow
D=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
L=/home/singjc/.jcode/scratch/dia2/ctest_cache/c/hela_entrap_codon.predicted_decoys.pqp
run() { n=$1; shift; /usr/bin/time -f "wall %E rss %M KB" $B -in $D -lib $L -out_dir $W/$n -threads 8 -seed 0 "$@" > $W/$n.log 2>&1; echo "JCODE_CHECKPOINT {\"message\":\"$n exit $? $(tail -1 $W/$n.log)\"}"; }
run ims06 -im_tolerance_scale 0.6 -features_out
run ims04 -im_tolerance_scale 0.4 -features_out
