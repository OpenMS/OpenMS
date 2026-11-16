#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH --mem=48G
#SBATCH -c 8
# PLAN_MS1: one self-calibrating OpenDIAWorkflow run, CUR (no flag) or NEW (-ms1_features).
# Usage: sbatch -J ms1-RR-sS-ARM -o logs/%x.%j.out run_ms1.sbatch RR S cur|new [extra flags]
# HOME is node-local: the Thermo reader's .NET named mutex deadlocked concurrent jobs on the shared NFS home.
C=/nfs/research/juan/singj/opendia
W=$C/rtwin2-20261116
R=$1; S=$2; A=$3; shift 3; O=$W/ms1/raw${R}_s${S}_$A
F=""; [ "$A" = new ] && F="-ms1_features"
mkdir -p $O
SC=${TMPDIR:-/tmp}/ms1-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $W/container/opendia-ms1-51fb7cc7ec.sif \
  OpenDIAWorkflow -in $C/proteobench-astral-entrapment/transfer/raw/LFQ_Astral_DIA_15min_50ng_Human_$R.raw \
  -lib $C/gap-20261114/home/lib/library.pqp -scoring_lib $W/lib/library12_full.pqp \
  -out_dir $O -threads 8 -seed $S $F "$@" > $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log; exit $rc
