#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH --mem=48G
#SBATCH -c 8
# PLAN_RTWIN2: one self-calibrating OpenDIAWorkflow run with -rt_width_selection (both arms from one run).
# Usage: sbatch -J rtwin-RR-sS -o logs/%x.%j.out run.sbatch RR S
C=/nfs/research/juan/singj/opendia
W=$C/rtwin2-20261116
R=$1; S=$2; O=$W/out/raw${R}_s$S
mkdir -p $O
SC=${TMPDIR:-/tmp}/rtwin-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v singularity exec --cleanenv --bind $C --bind $SC:/tmp $W/container/opendia-rtwin-d191f44f08.sif \
  OpenDIAWorkflow -in $C/proteobench-astral-entrapment/transfer/raw/LFQ_Astral_DIA_15min_50ng_Human_$R.raw \
  -lib $C/gap-20261114/home/lib/library.pqp -scoring_lib $W/lib/library12_full.pqp \
  -out_dir $O -threads 8 -seed $S -rt_width_selection > $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log; exit $rc
