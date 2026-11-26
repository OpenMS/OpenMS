#!/usr/bin/env bash
#SBATCH -J flscreen
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 08:00:00
#SBATCH -a 0-7
#SBATCH -o /nfs/research/juan/singj/opendia/foldlocal-20261126/logs/%A_%a.out
# Dev screen SCREEN_FOLDLOCAL.md: arms P (pooled) and L (fold_local) on calsearch D/H G0f seeds 0,1 (inputs read-only).
set -euo pipefail
C=/nfs/research/juan/singj/opendia; W=$C/foldlocal-20261126
RUNS=(D_G0f_s0 D_G0f_s1 H_G0f_s0 H_G0f_s1)
R=${RUNS[$((SLURM_ARRAY_TASK_ID / 2))]}; A=$([ $((SLURM_ARRAY_TASK_ID % 2)) = 0 ] && echo P || echo L)
S=${R##*_s}; X=; [ $A = L ] && X=fold_local
F=$(ls $C/calsearch-20261126/runs/$R/*.features.tsv)
/usr/bin/time -v $W/gbm_replay_static $F $W/hela_keys_e.tsv $S 8 $W/out/${R}_$A.tsv $X 2> $W/out/${R}_$A.log
