#!/bin/bash
#SBATCH --job-name=art2-evaluate
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -euo pipefail
export PATH=/ebi/slurm/codon/bin:$PATH
WORK=/nfs/research/juan/singj/opendia/astral-rt2-20261125
B=/nfs/research/juan/singj/opendia/astral-full-20261123
C=/nfs/research/juan/singj/opendia/astral-parquet-20261124
cd "$WORK"
mkdir -p results
sacct -j "$(cat primary_job_id.txt)" --parsable2 --units=K --format=JobID,State,ExitCode,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS,NodeList > results/slurm_accounting.tsv
set +e
python3 compare_resources.py "$WORK" "$B" > results/comparison.json
primary_rc=$?
python3 compare_resources.py "$WORK" "$C" > results/comparison_compact.json
compact_rc=$?
set -e
printf 'primary=%s\ncompact_reference=%s\n' "$primary_rc" "$compact_rc" > results/evaluation_exit_codes.txt
test "$primary_rc" = 0
