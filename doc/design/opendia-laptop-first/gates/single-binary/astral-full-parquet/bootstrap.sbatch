#!/bin/bash
#SBATCH --job-name=aparquet-bootstrap
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=24:00:00
set -euo pipefail
WORK=/nfs/research/juan/singj/opendia/astral-parquet-20261124
B=/nfs/research/juan/singj/opendia/astral-full-20261123
CONTAINER=/nfs/research/juan/singj/opendia/opendia-compact-parquet-20261124.sif
cd "$WORK"
test ! -e bootstrap
test ! -e library/search6.opendia.parquet
test ! -e library/scoring12.opendia.parquet
sha256sum -c preflight.sha256
mkdir bootstrap
scratch=$(mktemp -d "${TMPDIR:-/tmp}/aparquet-bootstrap-${SLURM_JOB_ID}-XXXXXX")
trap 'rm -rf "$scratch"' EXIT
{ hostname; lscpu; } > bootstrap/node.txt
set +e
/usr/bin/time -v -o bootstrap/time.txt singularity exec --cleanenv --no-home --home "$scratch:/tmp" --env TMPDIR=/tmp --bind /nfs/research/juan/singj/opendia --bind "$scratch:/tmp" "$CONTAINER" OpenDIA \
 -in /nfs/research/juan/singj/opendia/proteobench-astral-entrapment/transfer/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw \
 -lib "$WORK/library/search6.pqp" -scoring_lib "$WORK/library/scoring12.pqp" \
 -out_dir "$WORK/bootstrap" -threads 16 -seed 0 > bootstrap/run.log 2>&1
rc=$?
set -e
printf '%s\n' "$rc" > bootstrap/exit_code.txt
test "$rc" = 0
for name in precursors.tsv LFQ_Astral_DIA_15min_50ng_Human_01.precursors.tsv LFQ_Astral_DIA_15min_50ng_Human_01.calibration.json; do
 cmp "bootstrap/$name" "$B/runs/raw01_s0/$name"
done
test "$(cat "$B/runs/raw01_s0/exit_code.txt")" = 0
test -s library/search6.opendia.parquet
test -s library/scoring12.opendia.parquet
sha256sum library/search6.opendia.parquet library/scoring12.opendia.parquet > bootstrap/library.sha256
