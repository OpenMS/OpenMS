#!/bin/bash
W=/home/singjc/.jcode/scratch/dia2/gap7/ev
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIAWorkflow
D=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
L=/home/singjc/.jcode/scratch/dia2/ctest_cache/c/hela_entrap_codon.predicted_decoys.pqp
run() { n=$1; shift; /usr/bin/time -f "wall %E rss %M KB" $B -in $D -lib $L -out_dir $W/$n -threads 8 -seed 0 "$@" > $W/$n.log 2>&1; echo "JCODE_CHECKPOINT {\"message\":\"$n exit $? $(tail -1 $W/$n.log)\"}"; }
run ims10 -im_tolerance_scale 1.0
f=20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.precursors.tsv
cmp $W/ims10/$f $W/def2/$f && echo "ims10 SAME as default"
run ims08 -im_tolerance_scale 0.8
run ims05 -im_tolerance_scale 0.5
