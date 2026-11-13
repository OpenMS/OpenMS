#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --mem=96G
#SBATCH -c 12
# PLAN_DECOYGEN data job: the PAIRLIB chain_raw.sh, unchanged, on arm <X>'s run root (bound onto the PAIRLIB path).
# Usage: sbatch -J dg-<X>-<NN> -o logs/raw.%x.%j.out run_raw.sbatch <X> 01|02
DG=/nfs/research/juan/singj/opendia/decoygen-20261113
W=/nfs/research/juan/singj/opendia/pairlib-20261113
T=/nfs/research/juan/singj/opendia/proteobench-astral-entrapment/transfer
LOC=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113
RAWLOC=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw
SC=${TMPDIR:-/tmp}/dg-$SLURM_JOB_ID; mkdir -p $SC
singularity exec --cleanenv --bind $DG/arm$1/home:$LOC --bind $T/raw:$RAWLOC --bind $SC:/tmp $W/container/opendia-pairlib-a8be5bed3c.sif \
  bash $LOC/chain_raw.sh $2
rc=$?; rm -rf $SC; exit $rc
