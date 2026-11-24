#!/bin/bash
#SBATCH --job-name=aparquet-run
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --array=0-14%3
set -euo pipefail
export PATH=/ebi/slurm/codon/bin:$PATH
WORK=/nfs/research/juan/singj/opendia/astral-parquet-20261124
B=/nfs/research/juan/singj/opendia/astral-full-20261123
CONTAINER=/nfs/research/juan/singj/opendia/opendia-compact-parquet-20261124.sif
cd "$WORK"
sha256sum -c runtime.sha256
sha256sum -c bootstrap/library.sha256
# afterok plus actual process exit and byte identity, not a sentinel.
# Slurm accounting can lag job completion. It remains mandatory at final evaluation.
test "$(cat bootstrap/exit_code.txt)" = 0
for name in precursors.tsv LFQ_Astral_DIA_15min_50ng_Human_01.precursors.tsv LFQ_Astral_DIA_15min_50ng_Human_01.calibration.json; do
 cmp "bootstrap/$name" "$B/runs/raw01_s0/$name"
done
raw_index=$((SLURM_ARRAY_TASK_ID / 5 + 1))
seed=$((SLURM_ARRAY_TASK_ID % 5))
printf -v raw '%02d' "$raw_index"
out="$WORK/runs/raw${raw}_s${seed}"
test ! -e "$out"
mkdir -p "$out"
scratch=$(mktemp -d "${TMPDIR:-/tmp}/aparquet-${SLURM_JOB_ID}-XXXXXX")
trap 'rm -rf "$scratch"' EXIT
{ hostname; lscpu; printf 'raw=%s seed=%s threads=16 memory_budget=256G\n' "$raw" "$seed"; } > "$out/node.txt"
set +e
/usr/bin/time -v -o "$out/time.txt" singularity exec --cleanenv --no-home --home "$scratch:/tmp" --env TMPDIR=/tmp --bind /nfs/research/juan/singj/opendia --bind "$scratch:/tmp" "$CONTAINER" OpenDIA \
 -in "/nfs/research/juan/singj/opendia/proteobench-astral-entrapment/transfer/raw/LFQ_Astral_DIA_15min_50ng_Human_${raw}.raw" \
 -lib "$WORK/library/search6.opendia.parquet" -scoring_lib "$WORK/library/scoring12.opendia.parquet" \
 -out_dir "$out" -threads 16 -seed "$seed" > "$out/run.log" 2>&1
rc=$?
set -e
printf '%s\n' "$rc" > "$out/exit_code.txt"
exit "$rc"
