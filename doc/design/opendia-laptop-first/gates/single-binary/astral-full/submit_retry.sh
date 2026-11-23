#!/bin/bash
# LABEL_AMENDMENT.md attempt 3: views, then the unchanged 15-run array and evaluation. Run once.
set -euo pipefail
export PATH=/ebi/slurm/codon/bin:$PATH
cd /nfs/research/juan/singj/opendia/astral-full-20261123
test ! -e attempt2
test ! -e runs
sha256sum -c PRECHECK_ASTRAL_FULL.sha256
sha256sum -c diagnostics/attempt2_library.sha256
# Retain attempt 2: logs, accounting, job IDs and the failed view outputs.
mkdir attempt2
sacct -j "$(paste -sd, job_ids.txt)" --parsable2 --units=K --format=JobID,JobName,State,ExitCode,Submit,Start,End,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS,NodeList > attempt2/accounting.tsv
mv job_ids.txt submission_attempt2.log results attempt2/
cp -a logs attempt2/logs
mkdir attempt2/library
mv library/structure.json library/views.time.txt library/prepare.help.txt attempt2/library/
cp -a library/decoys.log library/decoys.time.txt attempt2/library/
views=$(sbatch --parsable --output=logs/views-%j.log prepare_views.sbatch)
printf '%s\n' "$views" > job_ids.txt
search=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$views" --output=logs/search-%A_%a.log run.sbatch)
printf '%s\n' "$search" >> job_ids.txt
evaluation=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterany:$search" --output=logs/evaluate-%j.log evaluate.sbatch)
printf '%s\n' "$evaluation" >> job_ids.txt
printf 'views=%s search=%s evaluation=%s\n' "$views" "$search" "$evaluation" | tee submission_attempt3.log
