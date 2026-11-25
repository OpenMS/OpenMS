#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=160G
#SBATCH -t 8:00:00
# PRECHECK_BLEND HeLa cell. Usage: sbatch -J bh-SET-ARM-sS -o logs/%x.%j.out blend_hela.sbatch D|H S ARM   (ARM: B5 | G0)
C=/nfs/research/juan/singj/opendia; W=$C/blend-20261125; V=$C/entlib-hela-20261118/data
DD=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
HD=$V/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
LIB=$C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp
P=$1; S=$2; A=$3
case $P in D) IN=$DD;; H) IN=$HD;; *) echo "unknown set $P"; exit 2;; esac
case $A in
  G0) F="";;
  B5) F="-gbm_blend 0.5";;
  *) echo "unknown arm $A"; exit 2;;
esac
O=$W/runs/${P}_${A}_s$S; test ! -e $O; mkdir -p $O
SC=${TMPDIR:-/tmp}/bh-$SLURM_JOB_ID; mkdir -p $SC
echo "node $(hostname) cpu $(lscpu | grep 'Model name' | sed 's/.*: *//') flags $F" > $O.log
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $C/opendia-blend-20261125.sif OpenDIA -in $IN -lib $LIB -out_dir $O -threads 8 -seed $S $F >> $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log
REF=$C/entlib-gbt-20261122/runs/${P}_gbt_s$S
if [ "$A" = G0 ] && [ -d $REF ]; then for f in $(cd $REF && ls); do cmp -s $O/$f $REF/$f && echo "GATE SAME $f" >> $O.log || echo "GATE DIFF $f" >> $O.log; done; fi
exit $rc
