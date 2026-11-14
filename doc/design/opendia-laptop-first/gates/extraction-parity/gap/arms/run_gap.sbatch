#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH --mem=96G
#SBATCH -c 12
# PLAN_GAP job: chain_gap.sh inside the PAIRLIB SIF on the gap run root. Usage: sbatch -J gap -o logs/gap.%j.out run_gap.sbatch
G=/nfs/research/juan/singj/opendia/gap-20261114
W=/nfs/research/juan/singj/opendia/pairlib-20261113
LOC=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113
PILOT=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110
REF=/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332
singularity exec --cleanenv --bind $G/home:$LOC --bind $W/ref:$REF --bind $W/pilot:$PILOT $W/container/opendia-pairlib-a8be5bed3c.sif \
  bash $LOC/chain_gap.sh
