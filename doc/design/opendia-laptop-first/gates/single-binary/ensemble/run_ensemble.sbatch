#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH --mem=48G
#SBATCH -c 8
# PLAN_ENSEMBLE: one self-calibrating OpenDIAWorkflow run with -ensemble M (CUR = 1, NEW = 3).
# Usage: sbatch -J ens-RR-sS-mM -o logs/%x.%j.out run_ensemble.sbatch RR S M
# HOME is the node-local job scratch: the Thermo reader's .NET named mutex lives under HOME, and on the shared NFS home
# concurrent jobs reading the same raw file deadlocked (RESULT_RTPROFILE infrastructure note).
C=/nfs/research/juan/singj/opendia
W=$C/rtwin2-20261116
R=$1; S=$2; M=$3; O=$W/ensemble/raw${R}_s${S}_m$M
mkdir -p $O
SC=${TMPDIR:-/tmp}/ens-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $W/container/opendia-ens-09aca5207c.sif \
  OpenDIAWorkflow -in $C/proteobench-astral-entrapment/transfer/raw/LFQ_Astral_DIA_15min_50ng_Human_$R.raw \
  -lib $C/gap-20261114/home/lib/library.pqp -scoring_lib $W/lib/library12_full.pqp \
  -out_dir $O -threads 8 -seed $S -ensemble $M > $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log; exit $rc
