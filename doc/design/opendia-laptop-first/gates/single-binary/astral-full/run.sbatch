#!/bin/bash
#SBATCH --job-name=afull-run
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --array=0-14%3
set -euo pipefail
cd /nfs/research/juan/singj/opendia/astral-full-20261123
test -s library/preparation.complete
raw_index=$((SLURM_ARRAY_TASK_ID / 5 + 1))
seed=$((SLURM_ARRAY_TASK_ID % 5))
printf -v raw '%02d' "$raw_index"
out="/nfs/research/juan/singj/opendia/astral-full-20261123/runs/raw${raw}_s${seed}"
test ! -e "$out"
mkdir -p "$out"
scratch=$(mktemp -d "${TMPDIR:-/tmp}/afull-run-${SLURM_JOB_ID}-XXXXXX")
{ hostname; lscpu; printf 'raw=%s seed=%s threads=16 memory_budget=256G\n' "$raw" "$seed"; } > "$out/node.txt"
set +e
/usr/bin/time -v -o "$out/time.txt" singularity exec --cleanenv --no-home --home "$scratch:/tmp" --env TMPDIR=/tmp --bind /nfs/research/juan/singj/opendia --bind "$scratch:/tmp" /nfs/research/juan/singj/opendia/opendia-pd10-fastint2.sif OpenDIA \
 -in "/nfs/research/juan/singj/opendia/proteobench-astral-entrapment/transfer/raw/LFQ_Astral_DIA_15min_50ng_Human_${raw}.raw" \
 -lib /nfs/research/juan/singj/opendia/astral-full-20261123/library/search6.pqp \
 -scoring_lib /nfs/research/juan/singj/opendia/astral-full-20261123/library/scoring12.pqp \
 -out_dir "$out" -threads 16 -seed "$seed" > "$out/run.log" 2>&1
rc=$?
set -e
printf '%s\n' "$rc" > "$out/exit_code.txt"
exit "$rc"
