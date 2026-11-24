#!/bin/bash
#SBATCH --job-name=aparquet-evaluate
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -euo pipefail
export PATH=/ebi/slurm/codon/bin:$PATH
WORK=/nfs/research/juan/singj/opendia/astral-parquet-20261124
B=/nfs/research/juan/singj/opendia/astral-full-20261123
cd "$WORK"
mkdir -p results
sacct -j "$(cat bootstrap_job_id.txt)" --parsable2 --units=K --format=JobID,State,ExitCode,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS,NodeList > bootstrap_accounting.tsv
sacct -j "$(cat primary_job_id.txt)" --parsable2 --units=K --format=JobID,State,ExitCode,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS,NodeList > results/slurm_accounting.tsv
set +e
python3 compare_resources.py "$WORK" "$B" --bootstrap > results/bootstrap.json
bootstrap_rc=$?
python3 compare_resources.py "$WORK" "$B" > results/comparison.json
primary_rc=$?
set -e
printf 'bootstrap=%s\nprimary=%s\n' "$bootstrap_rc" "$primary_rc" > results/evaluation_exit_codes.txt
test "$bootstrap_rc" = 0 && test "$primary_rc" = 0
