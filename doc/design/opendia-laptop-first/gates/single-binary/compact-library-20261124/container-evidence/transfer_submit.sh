#!/usr/bin/env bash
set -euo pipefail
repo=/home/singjc/Documents/github/OpenMS-wp00-design-docs
bundle="$repo/doc/design/opendia-laptop-first/gates/single-binary/astral-full-parquet"
work=/home/singjc/.jcode/scratch/compact-20261124
sif="$work/container/opendia-compact-parquet-20261124.sif"
ssh_options=(-o BatchMode=yes -o ConnectTimeout=20 -o ControlMaster=no -o ControlPath=none -o ServerAliveInterval=30 -o ServerAliveCountMax=3)
# Refuse uncommitted or unpushed approval. Existing remote files are never replaced.
git -C "$repo" cat-file -e HEAD:doc/design/opendia-laptop-first/gates/single-binary/astral-full-parquet/PREFLIGHT_APPROVED.md
git -C "$repo" diff --quiet HEAD -- "$bundle"
head=$(git -C "$repo" rev-parse HEAD)
pushed=$(GIT_SSH_COMMAND='ssh -o BatchMode=yes -o ConnectTimeout=20 -o ControlMaster=no -o ControlPath=none' git -C "$repo" ls-remote origin refs/heads/docs/opendia-wp00-design-plan | cut -f1)
test "$head" = "$pushed"
printf 'Pushed approval commit: %s\n' "$head"
sha256sum "$sif"
# One unique partial file is retained on transport failure, not removed or retried blindly.
timeout 1800 ssh "${ssh_options[@]}" codon "become pst_prd bash -c 'set -euo pipefail; target=/nfs/research/juan/singj/opendia/opendia-compact-parquet-20261124.sif; test ! -e \"\$target\"; test ! -L \"\$target\"; test ! -e \"\$target.partial\"; test ! -L \"\$target.partial\"; umask 002; cat > \"\$target.partial\"; printf \"035cd7004896d636dc83800230876b4a504fbc247433ea180c02e7ddaffb9caf  %s.partial\\n\" \"\$target\" | sha256sum -c -; mv -T \"\$target.partial\" \"\$target\"'" < "$sif"
printf 'JCODE_PROGRESS {"message":"unique Codon SIF transferred and hash verified"}\n'
# Tar is locally constructed from the eleven exact approved files, no result directories.
tar -cf - -C "$bundle" PLAN_PARQUET.md PREFLIGHT_APPROVED.md source_provenance.txt preflight.sha256 runtime.sha256 bootstrap.sbatch run.sbatch evaluate.sbatch submit.sh compare_resources.py test_compare_resources.py |
 timeout 120 ssh "${ssh_options[@]}" codon "become pst_prd bash -c 'set -euo pipefail; work=/nfs/research/juan/singj/opendia/astral-parquet-20261124; test ! -e \"\$work\"; test ! -L \"\$work\"; umask 002; mkdir \"\$work\"; tar -xf - -C \"\$work\"'"
printf 'JCODE_PROGRESS {"message":"frozen bundle transferred into fresh Codon directory"}\n'
# The frozen submit script verifies full preflight before the first sbatch. Keep all output.
timeout 1800 ssh "${ssh_options[@]}" codon "become pst_prd bash -lc 'cd /nfs/research/juan/singj/opendia/astral-parquet-20261124 && bash submit.sh'"
