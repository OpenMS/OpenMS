#!/bin/bash
# PLAN_RTCENTRE design runs: HeLa D, seeds 0-4, BASE (0e9e4fc3a5) and NEW (6dd893b954); then Astral raw01 seed 0 NEW.
# Two runs at a time, 4 threads each. Output is thread-count independent (checked earlier), so 4 vs 8 threads is fine.
B=/home/singjc/.jcode/scratch/dia2
H=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published
D=$H/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
LIB=$H/200623_hela_evosep_fractions_library_ptypic_decoy.pqp
mkdir -p $B/rtc
run() { # $1 arm, $2 seed
  /usr/bin/time -v $B/bin/OpenDIAWorkflow.$1 -in $D -lib $LIB -out_dir $B/rtc/D_$1_s$2 -threads 4 -seed $2 > $B/rtc/D_$1_s$2.log 2>&1
  echo "D $1 s$2 $(grep -o 'q <= 0.01: [0-9]*' $B/rtc/D_$1_s$2.log)"
}
for s in 0 1 2 3 4; do
  run base $s & run rtcentre $s & wait
done
G=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
DOTNET_ROOT=/home/singjc/.dotnet/runtime8 $B/bin/OpenDIAWorkflow.rtcentre -in /home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw \
  -lib $G/gap-20261114/armE_raw01/library.pqp -scoring_lib $G/lib12full/pairlib-view/lib12/library12_full.pqp \
  -out_dir $B/rtc/A_rtcentre -threads 8 -seed 0 > $B/rtc/A_rtcentre.log 2>&1
DOTNET_ROOT=/home/singjc/.dotnet/runtime8 $B/bin/OpenDIAWorkflow.base -in /home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw \
  -lib $G/gap-20261114/armE_raw01/library.pqp -scoring_lib $G/lib12full/pairlib-view/lib12/library12_full.pqp \
  -out_dir $B/rtc/A_base -threads 8 -seed 0 > $B/rtc/A_base.log 2>&1
for f in precursors.tsv LFQ_Astral_DIA_15min_50ng_Human_01.precursors.tsv LFQ_Astral_DIA_15min_50ng_Human_01.calibration.json; do
  cmp $B/rtc/A_base/$f $B/rtc/A_rtcentre/$f && echo "A $f identical"
done
echo done
