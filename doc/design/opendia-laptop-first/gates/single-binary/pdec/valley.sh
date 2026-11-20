#!/bin/bash
# PRECHECK_VALLEY (diagnostic, seed 0 only): r~1.56 HeLa D, current default vs search picker settings.
# def: byte-identity of the new options at their defaults (compare with ev/def2, same default build output).
W=/home/singjc/.jcode/scratch/dia2/gap7/ev
B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIAWorkflow
D=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
L=/home/singjc/.jcode/scratch/dia2/ctest_cache/c/hela_entrap_codon.predicted_decoys.pqp
run() { n=$1; shift; /usr/bin/time -f "wall %E rss %M KB" $B -in $D -lib $L -out_dir $W/$n -threads 8 -seed 0 "$@" > $W/$n.log 2>&1; echo "JCODE_CHECKPOINT {\"message\":\"$n exit $? $(tail -1 $W/$n.log)\"}"; }
run vdef
f=20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.precursors.tsv
cmp $W/vdef/$f $W/def2/$f && echo "vdef SAME as def2"
run v08m8 -valley_fraction 0.8 -max_peaks 8 -features_out
run v05m8 -max_peaks 8
