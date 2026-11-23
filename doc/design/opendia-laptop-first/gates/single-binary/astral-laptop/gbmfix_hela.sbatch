#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=160G
#SBATCH -t 8:00:00
# PRECHECK_GBMFIX HeLa arm. Usage: sbatch -J gfh-ARM-sS -o logs/%x.%j.out gbmfix_hela.sbatch D7|D156 S ARM
# D7: r~7 D (hela_r5 library, as PLAN_GBT); D156: the r~1.56 D byte-check gate (G0 only).
C=/nfs/research/juan/singj/opendia; W=$C/gbmfix2-20261123
V=$C/entlib-hela-20261118/data
DD=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d; N=20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737
L=$1; S=$2; A=$3
case $A in
  G0) F="";;
  G1) F="-gbm_rounds 1";;
  GQ) F="-gbm_positive_q 0.005";;
  *) echo "unknown arm $A"; exit 2;;
esac
case $L in
  D7) LIB=$C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp; REF=$C/entlib-gbt-20261122/runs/D_gbt_s$S;;
  D156) LIB=/homes/sing/hela_entrap_codon.predicted_decoys.pqp; REF=$C/gbmdef-20261122/D_def_s0;;
esac
O=$W/runs/${L}_${A}_s$S; mkdir -p $O
SC=${TMPDIR:-/tmp}/gfh-$SLURM_JOB_ID; mkdir -p $SC
echo "node $(hostname) cpu $(lscpu | grep 'Model name' | sed 's/.*: *//') flags $F" > $O.log
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind /homes/sing --bind $SC:/tmp \
  $C/opendia-pd10-gbmvar.sif OpenDIA -in $DD -lib $LIB -out_dir $O -threads 8 -seed $S $F >> $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log
if [ "$A" = G0 ]; then for f in $(cd $REF && ls); do cmp -s $O/$f $REF/$f && echo "GATE SAME $f" >> $O.log || echo "GATE DIFF $f" >> $O.log; done; fi
exit $rc
