#!/usr/bin/env bash
#SBATCH -J bnscreen
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH -a 0-3
#SBATCH -o /nfs/research/juan/singj/opendia/bulknull-20261126/logs/%A_%a.out
# SCREEN_BULKNULL: bulkNullQ on the exported final score of calsearch D/H G0f seeds 0,1 (inputs read-only).
set -euo pipefail
C=/nfs/research/juan/singj/opendia; W=$C/bulknull-20261126
RUNS=(D_G0f_s0 D_G0f_s1 H_G0f_s0 H_G0f_s1)
R=${RUNS[$SLURM_ARRAY_TASK_ID]}
F=$(ls $C/calsearch-20261126/runs/$R/*.features.tsv)
/usr/bin/time -v $W/bulknull_replay $F $W/hela_keys_relabel.tsv 8 $W/out/$R.tsv 2> $W/out/$R.log
