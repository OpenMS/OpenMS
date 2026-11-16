#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH --mem=48G
#SBATCH -c 8
# PLAN_RTPROFILE: one self-calibrating OpenDIAWorkflow run with -id_recalibration (CUR = pass1, NEW = pass 2).
# Usage: sbatch -J rtp-RR-sS -o logs/%x.%j.out run_rtprofile.sbatch RR S
# HOME and TMPDIR point to the node-local job scratch. The Thermo reader's .NET named mutex (CBaseIO_<raw path>) lives
# under HOME; on the NFS home it is shared by all nodes, and concurrent jobs reading the same raw file deadlocked on it.
C=/nfs/research/juan/singj/opendia
W=$C/rtwin2-20261116
R=$1; S=$2; O=$W/rtprofile/raw${R}_s$S
mkdir -p $O
SC=${TMPDIR:-/tmp}/rtp-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp $W/container/opendia-rtprofile-f86fb0db22.sif \
  OpenDIAWorkflow -in $C/proteobench-astral-entrapment/transfer/raw/LFQ_Astral_DIA_15min_50ng_Human_$R.raw \
  -lib $C/gap-20261114/home/lib/library.pqp -scoring_lib $W/lib/library12_full.pqp \
  -out_dir $O -threads 8 -seed $S -id_recalibration > $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log; exit $rc
