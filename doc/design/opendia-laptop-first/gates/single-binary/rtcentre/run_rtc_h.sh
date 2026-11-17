#!/bin/bash
# PLAN_RTCENTRE held-out runs: HeLa H (22474.d), seeds 0-4, BASE and NEW, two at a time, 4 threads each.
B=/home/singjc/.jcode/scratch/dia2
H=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published
D=/home/singjc/Documents/github/pyopenms-viz/nbs/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
LIB=$H/200623_hela_evosep_fractions_library_ptypic_decoy.pqp
mkdir -p $B/rtc
run() {
  /usr/bin/time -v $B/bin/OpenDIAWorkflow.$1 -in $D -lib $LIB -out_dir $B/rtc/H_$1_s$2 -threads 4 -seed $2 > $B/rtc/H_$1_s$2.log 2>&1
  echo "H $1 s$2 $(grep -o 'q <= 0.01: [0-9]*' $B/rtc/H_$1_s$2.log)"
}
for s in 0 1 2 3 4; do
  run base $s & run rtcentre $s & wait
done
echo done
