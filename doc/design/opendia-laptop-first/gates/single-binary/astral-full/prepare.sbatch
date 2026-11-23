#!/bin/bash
#SBATCH --job-name=afull-prepare
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=24:00:00
set -euo pipefail
cd /nfs/research/juan/singj/opendia/astral-full-20261123
scratch=$(mktemp -d "${TMPDIR:-/tmp}/afull-prepare-${SLURM_JOB_ID}-XXXXXX")
test -s library/generation.complete
/usr/bin/time -v -o library/decoys.time.txt singularity exec --cleanenv --no-home --home "$scratch:/tmp" --bind /nfs/research/juan/singj/opendia --bind "$scratch:/tmp" /nfs/research/juan/singj/opendia/opendia-pd10-fastint2.sif OpenSwathDecoyGenerator \
 -in /nfs/research/juan/singj/opendia/astral-full-20261123/library/targets20.pqp \
 -out /nfs/research/juan/singj/opendia/astral-full-20261123/library/predicted20.pqp \
 -method pseudo-reverse -predict_decoys true -predict_keep_short true -predict_nce 25 -predict_instrument Lumos -threads 16 > library/decoys.log 2>&1
"$HOME/.local/bin/uv" run -q --no-project python prepare_library.py --help > library/prepare.help.txt
# Library view construction command is frozen with prepare_library.py before submission.
/usr/bin/time -v -o library/views.time.txt "$HOME/.local/bin/uv" run -q --no-project python prepare_library.py library/predicted20.pqp --search6 library/search6.pqp --scoring12 library/scoring12.pqp > library/structure.json
sha256sum library/predicted20.pqp library/search6.pqp library/scoring12.pqp > library/prepared.sha256
printf 'preparation_complete\n' > library/preparation.complete
