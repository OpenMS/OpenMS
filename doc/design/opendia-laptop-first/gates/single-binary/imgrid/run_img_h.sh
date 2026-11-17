#!/bin/bash
# PLAN_IMGRID held-out: H NEW seeds 0-4 (BASE reused from PLAN_RTCENTRE H_rtcentre_s*).
B=/home/singjc/.jcode/scratch/dia2
H=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published
D=/home/singjc/Documents/github/pyopenms-viz/nbs/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
run() {
  $B/bin/OpenDIAWorkflow.imgrid -in $D -lib $H/200623_hela_evosep_fractions_library_ptypic_decoy.pqp -out_dir $B/rtc/H_imgrid_s$1 -threads 4 -seed $1 > $B/rtc/H_imgrid_s$1.log 2>&1
  echo "H imgrid s$1 $(grep -o 'q <= 0.01: [0-9]*' $B/rtc/H_imgrid_s$1.log)"
}
run 0 & run 1 & wait
run 2 & run 3 & wait
run 4
echo done
