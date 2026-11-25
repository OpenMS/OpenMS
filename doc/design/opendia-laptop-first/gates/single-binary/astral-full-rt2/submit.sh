#!/bin/bash
# Execute only after the root's pushed preflight amendment, never during drafting.
set -euo pipefail
export PATH=/ebi/slurm/codon/bin:$PATH
WORK=/nfs/research/juan/singj/opendia/astral-rt2-20261125
C=/nfs/research/juan/singj/opendia/astral-parquet-20261124
cd "$WORK"
for path in job_ids.txt runs library primary_job_id.txt; do
 test ! -e "$path" && test ! -L "$path"
done
# Atomic reservation prevents simultaneous duplicate submissions. Never remove it automatically.
mkdir submission.lock
test -s PREFLIGHT_APPROVED.md
test -s preflight.sha256
test -s runtime.sha256
sha256sum -c preflight.sha256
sha256sum -c runtime.sha256
mkdir library logs
for name in search6 scoring12; do
 ln -s "$C/library/$name.opendia.parquet" "library/$name.opendia.parquet"
done
printf '%s\n' \
 '1943673979d97b8da9b87cd6b9df0a48ba5519d72af1537bf64f582f05dafaeb  library/search6.opendia.parquet' \
 '425c9e02a3e9fe0a21ccfce44d63b542e1a4664a6d5a9607b0e349a0310603fa  library/scoring12.opendia.parquet' > library.sha256
sha256sum -c library.sha256
primary=$(sbatch --parsable --output=logs/run-%A_%a.log run.sbatch)
primary=${primary%%;*}
printf '%s\n' "$primary" > job_ids.txt
printf '%s\n' "$primary" > primary_job_id.txt
evaluation=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterany:$primary" --output=logs/evaluate-%j.log evaluate.sbatch)
evaluation=${evaluation%%;*}
printf '%s\n' "$evaluation" >> job_ids.txt
printf 'primary=%s evaluation=%s\n' "$primary" "$evaluation"
