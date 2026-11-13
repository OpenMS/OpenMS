#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --mem=96G
#SBATCH -c 12
# Usage: sbatch -J pairlib-rawNN -o logs/raw%x.%j.out run_raw.sbatch 01|02
W=/nfs/research/juan/singj/opendia/pairlib-20261113
T=/nfs/research/juan/singj/opendia/proteobench-astral-entrapment/transfer
LOC=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113
RAWLOC=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw
SC=${TMPDIR:-/tmp}/pairlib-$SLURM_JOB_ID; mkdir -p $SC
singularity exec --cleanenv --bind $W/home:$LOC --bind $T/raw:$RAWLOC --bind $SC:/tmp $W/container/opendia-pairlib-a8be5bed3c.sif \
  bash $LOC/chain_raw.sh $1
rc=$?; rm -rf $SC; exit $rc
