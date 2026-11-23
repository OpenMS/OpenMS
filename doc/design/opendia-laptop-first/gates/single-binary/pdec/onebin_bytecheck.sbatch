#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -t 6:00:00
# Cleanup byte check. Usage: run.sbatch D156|H7|D7off
C=/nfs/research/juan/singj/opendia; W=$C/onebin-20261123; V=$C/entlib-hela-20261118/data
SIF=$C/opendia-pd10-onebin.sif; A=$1; F=""
DD=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d; HD=$V/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
R7=$C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp
case $A in
  D156) IN=$DD; LIB=/homes/sing/hela_entrap_codon.predicted_decoys.pqp; REF=$C/gbmdef-20261122/D_def_s0;;
  D7V) IN=$DD; LIB=$R7; REF=$C/v7m6t-20261122/D_r7_V7M6_1; F="-valley_fraction 0.7 -max_peaks 6";;
  H7) IN=$HD; LIB=$R7; REF=$C/entlib-gbt-20261122/runs/H_gbt_s0;;
esac
O=$W/$A; mkdir -p $O; SC=${TMPDIR:-/tmp}/fi-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind /homes/sing --bind $SC:/tmp \
  $SIF OpenDIAWorkflow -in $IN -lib $LIB -out_dir $O -threads 8 -seed 0 $F > $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log
for f in $(cd $REF && ls); do cmp -s $O/$f $REF/$f && echo "SAME $f" >> $O.log || echo "DIFF $f" >> $O.log; done
exit $rc
