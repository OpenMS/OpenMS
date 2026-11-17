#!/bin/bash
# PLAN_IMGRID design runs: D NEW seeds 0-4 (BASE reused from PLAN_RTCENTRE D_rtcentre_s*), Astral raw01 identity.
B=/home/singjc/.jcode/scratch/dia2
H=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published
D=$H/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
LIB=$H/200623_hela_evosep_fractions_library_ptypic_decoy.pqp
run() {
  $B/bin/OpenDIAWorkflow.imgrid -in $D -lib $LIB -out_dir $B/rtc/D_imgrid_s$1 -threads 4 -seed $1 > $B/rtc/D_imgrid_s$1.log 2>&1
  echo "D imgrid s$1 $(grep -o 'q <= 0.01: [0-9]*' $B/rtc/D_imgrid_s$1.log)"
}
run 0 & run 1 & wait
run 2 & run 3 & wait
run 4
G=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
DOTNET_ROOT=/home/singjc/.dotnet/runtime8 $B/bin/OpenDIAWorkflow.imgrid -in /home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw \
  -lib $G/gap-20261114/armE_raw01/library.pqp -scoring_lib $G/lib12full/pairlib-view/lib12/library12_full.pqp \
  -out_dir $B/rtc/A_imgrid -threads 8 -seed 0 > $B/rtc/A_imgrid.log 2>&1
for f in precursors.tsv LFQ_Astral_DIA_15min_50ng_Human_01.precursors.tsv LFQ_Astral_DIA_15min_50ng_Human_01.calibration.json; do
  cmp $B/rtc/A_rtcentre/$f $B/rtc/A_imgrid/$f && echo "A $f identical"
done
echo done
