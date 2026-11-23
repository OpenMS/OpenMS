#!/bin/bash
#SBATCH --job-name=afull-views
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=24:00:00
# LABEL_AMENDMENT.md: only the view step of prepare.sbatch, on the unmodified attempt-2 libraries.
set -euo pipefail
cd /nfs/research/juan/singj/opendia/astral-full-20261123
test -s library/generation.complete
sha256sum -c diagnostics/attempt2_library.sha256
test ! -e library/search6.pqp
test ! -e library/scoring12.pqp
"$HOME/.local/bin/uv" run -q --no-project python prepare_library.py --help > library/prepare.help.txt
/usr/bin/time -v -o library/views.time.txt "$HOME/.local/bin/uv" run -q --no-project python prepare_library.py library/predicted20.pqp --search6 library/search6.pqp --scoring12 library/scoring12.pqp > library/structure.json
sha256sum library/predicted20.pqp library/search6.pqp library/scoring12.pqp > library/prepared.sha256
printf 'preparation_complete\n' > library/preparation.complete
