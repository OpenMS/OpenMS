#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 4:00:00
# PRECHECK_ASTQ ladder arm. Usage: sbatch -J astq-RR-ARM-sS -o logs/%x.%j.out astq.sbatch RR S ARM
# Container opendia-pd10-gbmdef.sif (OpenMS bcb156e036, the GBM default before the cleanup; still has the stage switches).
C=/nfs/research/juan/singj/opendia; W=$C/astq-20261123
R=$1; S=$2; A=$3
case $A in
  A0) F="";;
  L1) F="-no_gbm_rescore";;
  L2) F="-no_gbm_rescore -no_sibling_features";;
  L3) F="-no_gbm_rescore -no_sibling_features -no_rt_refine";;
  L4) F="-no_gbm_rescore -no_sibling_features -no_rt_refine -no_matched_filter";;
  L5) F="-no_gbm_rescore -no_sibling_features -no_rt_refine -no_matched_filter -trace_tolerance_ppm 0 -im_tolerance_scale 1 -no_extra_scores";;
  *) echo "unknown arm $A"; exit 2;;
esac
O=$W/runs/raw${R}_${A}_s$S; mkdir -p $O
SC=${TMPDIR:-/tmp}/astq-$SLURM_JOB_ID; mkdir -p $SC
echo "node $(hostname) cpu $(lscpu | grep 'Model name' | sed 's/.*: *//') flags $F" > $O.log
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $C/opendia-pd10-gbmdef.sif \
  OpenDIAWorkflow -in $C/proteobench-astral-entrapment/transfer/raw/LFQ_Astral_DIA_15min_50ng_Human_$R.raw \
  -lib $C/gap-20261114/home/lib/library.pqp -scoring_lib $C/rtwin2-20261116/lib/library12_full.pqp \
  -out_dir $O -threads 16 -seed $S $F >> $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log; exit $rc
