#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH --mem=96G
#SBATCH -c 12
W=/nfs/research/juan/singj/opendia/pairlib-20261113
T=/nfs/research/juan/singj/opendia/proteobench-astral-entrapment/transfer
LOC=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113
REF=/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332
singularity exec --cleanenv --bind $W/home:$LOC --bind $W/ref:$REF $W/container/opendia-pairlib-a8be5bed3c.sif bash $LOC/chain_score.sh
