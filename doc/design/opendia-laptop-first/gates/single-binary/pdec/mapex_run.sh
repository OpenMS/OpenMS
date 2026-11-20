#!/bin/sh
# usage: run.sh <name> [extra args]   (main build, r~1.56 D, 8 threads, nice)
N=$1; shift
O=/home/singjc/.jcode/scratch/dia2/mapex/$N; rm -rf $O; mkdir -p $O
nice -n 10 /usr/bin/time -v /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIAWorkflow -in /home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d -lib /home/singjc/.jcode/scratch/dia2/ctest_cache/c/hela_entrap_codon.predicted_decoys.pqp -out_dir $O -threads 8 -seed 0 "$@" > $O/log.txt 2>&1
echo exit $? >> $O/log.txt
