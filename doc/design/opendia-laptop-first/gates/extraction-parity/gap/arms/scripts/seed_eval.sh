#!/bin/bash
#SBATCH -t 03:00:00
#SBATCH --mem=48G
#SBATCH -c 12
# Post-hoc seed check: offsets 1..4 for D11FD and D11FDI, then compare_gap per offset.
G=/nfs/research/juan/singj/opendia/gap-20261114
W=/nfs/research/juan/singj/opendia/pairlib-20261113
LOC=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113
PILOT=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110
REF=/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python
S="singularity exec --cleanenv --bind $G/home:$LOC --bind $W/ref:$REF --bind $W/pilot:$PILOT $W/container/opendia-pairlib-a8be5bed3c.sif"
for o in 1 2 3 4; do
  $S $PY $LOC/scripts_gap/seed_check.py $o > $G/home/score/seed$o.log 2>&1 &
done
wait
for o in 1 2 3 4; do
  $S $PY $LOC/scripts_gap/compare_gap.py D11FD_s${o}WT D11FDI_s${o}WT > $G/home/score/SEED_COMPARE_$o.txt 2>> $G/home/score/seed_eval.err
done
