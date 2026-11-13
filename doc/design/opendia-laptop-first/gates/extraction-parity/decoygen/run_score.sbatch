#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH --mem=96G
#SBATCH -c 12
# PLAN_DECOYGEN scoring job: the PAIRLIB chain_score.sh, unchanged, on arm <X>'s run root, then the label-free null rates
# (deployable pre-check), the labelled null rates and the pooled rho (both evaluation only).
# Usage: sbatch -J dg-<X>-score -d afterok:<raw jobs> -o logs/score.%x.%j.out run_score.sbatch <X>
DG=/nfs/research/juan/singj/opendia/decoygen-20261113
W=/nfs/research/juan/singj/opendia/pairlib-20261113
LOC=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113
PILOT=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110
REF=/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python
S="singularity exec --cleanenv --bind $DG/arm$1/home:$LOC --bind $W/ref:$REF --bind $W/pilot:$PILOT $W/container/opendia-pairlib-a8be5bed3c.sif"
$S bash $LOC/chain_score.sh; rc=$?
$S $PY $LOC/scripts/null_rates.py > $DG/arm$1/home/score/NULLRATES.txt 2>&1
$S $PY $LOC/scripts/null_rates.py --labels > $DG/arm$1/home/score/NULLRATES_labels.txt 2>&1
[ -s $DG/arm$1/home/score/DIAG_D11FDWT.txt ] && $S $PY $LOC/scripts/pooled_rho.py $LOC D11FDWT > $DG/arm$1/home/score/RHO_D11FDWT.txt 2>&1
[ -s $DG/arm$1/home/score/DIAG_D11FWT.txt ] && $S $PY $LOC/scripts/pooled_rho.py $LOC D11FWT > $DG/arm$1/home/score/RHO_D11FWT.txt 2>&1
exit $rc
