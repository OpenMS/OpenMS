#!/bin/bash
# Execute only after the root's pushed preflight amendment, never during drafting.
set -euo pipefail
export PATH=/ebi/slurm/codon/bin:$PATH
WORK=/nfs/research/juan/singj/opendia/astral-parquet-20261124
B=/nfs/research/juan/singj/opendia/astral-full-20261123
cd "$WORK"
for path in job_ids.txt runs bootstrap bootstrap_job_id.txt primary_job_id.txt; do
 test ! -e "$path" && test ! -L "$path"
done
# Atomic reservation prevents simultaneous duplicate submissions. Never remove it automatically.
mkdir submission.lock
test -s PREFLIGHT_APPROVED.md
test -s preflight.sha256
test -s runtime.sha256
sha256sum -c preflight.sha256
sha256sum -c runtime.sha256
mkdir -p library logs
for name in search6 scoring12; do
 test ! -e "library/$name.pqp" && test ! -L "library/$name.pqp"
 ln -s "$B/library/$name.pqp" "library/$name.pqp"
 test ! -e "library/$name.opendia.parquet"
done
printf '%s\n' \
 '200dea366ea358f9f97674715c070a6bd4ade8a839d2310fbf214fbf845847c7  library/search6.pqp' \
 '95a19dec2d169d6b19f973b27eda9763df26fcad69d1a3847115cb76186c49a0  library/scoring12.pqp' | sha256sum -c -
bootstrap=$(sbatch --parsable --output=logs/bootstrap-%j.log bootstrap.sbatch)
bootstrap=${bootstrap%%;*}
printf '%s\n' "$bootstrap" > job_ids.txt
printf '%s\n' "$bootstrap" > bootstrap_job_id.txt
primary=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$bootstrap" --output=logs/run-%A_%a.log run.sbatch)
primary=${primary%%;*}
printf '%s\n' "$primary" >> job_ids.txt
printf '%s\n' "$primary" > primary_job_id.txt
evaluation=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterany:$primary" --output=logs/evaluate-%j.log evaluate.sbatch)
evaluation=${evaluation%%;*}
printf '%s\n' "$evaluation" >> job_ids.txt
printf 'bootstrap=%s primary=%s evaluation=%s\n' "$bootstrap" "$primary" "$evaluation"
