#!/usr/bin/env bash
# Guarded transfer and submission for the frozen rt2 experiment. Same structure as compact-20261124/transfer_submit.sh.
set -euo pipefail
repo=/home/singjc/Documents/github/OpenMS-wp00-design-docs
bundle="$repo/doc/design/opendia-laptop-first/gates/single-binary/astral-full-rt2"
work=/home/singjc/.jcode/scratch/rt2-20261125
sif="$work/container/opendia-rt2-20261125.sif"
ssh_options=(-o BatchMode=yes -o ConnectTimeout=20 -o ControlMaster=no -o ControlPath=none -o ServerAliveInterval=30 -o ServerAliveCountMax=3)
# Refuse uncommitted or unpushed approval. Existing remote files are never replaced.
git -C "$repo" cat-file -e HEAD:doc/design/opendia-laptop-first/gates/single-binary/astral-full-rt2/PREFLIGHT_APPROVED.md
git -C "$repo" diff --quiet HEAD -- "$bundle"
head=$(git -C "$repo" rev-parse HEAD)
pushed=$(GIT_SSH_COMMAND='ssh -o BatchMode=yes -o ConnectTimeout=20 -o ControlMaster=no -o ControlPath=none' git -C "$repo" ls-remote origin refs/heads/docs/opendia-wp00-design-plan | cut -f1)
test "$head" = "$pushed"
printf 'Pushed approval commit: %s\n' "$head"
sha256sum "$sif"
# One unique partial file is retained on transport failure, not removed or retried blindly.
timeout 1800 ssh "${ssh_options[@]}" codon "become pst_prd bash -c 'set -euo pipefail; target=/nfs/research/juan/singj/opendia/opendia-rt2-20261125.sif; test ! -e \"\$target\"; test ! -L \"\$target\"; test ! -e \"\$target.partial\"; test ! -L \"\$target.partial\"; umask 002; cat > \"\$target.partial\"; printf \"775a7bcd493820b513abee1b4f81b2f7e692e5862649009e583fd0be17092e1f  %s.partial\\n\" \"\$target\" | sha256sum -c -; mv -T \"\$target.partial\" \"\$target\"'" < "$sif"
printf 'JCODE_PROGRESS {"message":"unique Codon SIF transferred and hash verified"}\n'
# Tar is locally constructed from the ten exact approved files, no result directories.
tar -cf - -C "$bundle" PLAN_RT2.md PREFLIGHT_APPROVED.md source_provenance.txt preflight.sha256 runtime.sha256 run.sbatch evaluate.sbatch submit.sh compare_resources.py test_compare_resources.py |
 timeout 120 ssh "${ssh_options[@]}" codon "become pst_prd bash -c 'set -euo pipefail; work=/nfs/research/juan/singj/opendia/astral-rt2-20261125; test ! -e \"\$work\"; test ! -L \"\$work\"; umask 002; mkdir \"\$work\"; tar -xf - -C \"\$work\"'"
printf 'JCODE_PROGRESS {"message":"frozen bundle transferred into fresh Codon directory"}\n'
# The frozen submit script verifies full preflight before the first sbatch. Keep all output.
timeout 1800 ssh "${ssh_options[@]}" codon "become pst_prd bash -lc 'cd /nfs/research/juan/singj/opendia/astral-rt2-20261125 && bash submit.sh'"
