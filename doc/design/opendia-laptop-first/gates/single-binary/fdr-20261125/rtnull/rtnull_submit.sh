#!/usr/bin/env bash
set -euo pipefail
export PATH=/ebi/slurm/codon/bin:$PATH
W=/nfs/research/juan/singj/opendia/rtnull-20261126
cd "$W/frozen"
sha256sum -c PRECHECK_RTNULL.sha256
printf '%s  %s\n' c43320a6e0d350efa29d80db2cafe2470bd026edd5a4fbafe6eba6f5792f191a /nfs/research/juan/singj/opendia/evshare-20261126/opendia-evshare-20261126.sif | sha256sum -c -
printf '%s  %s\n' 76f2a6c22277eaf464a301b446e97b0fd909fbf6b26d79736b1b0337c304f15c /nfs/research/juan/singj/opendia/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp | sha256sum -c -
test ! -e "$W/submitted.tsv"
mkdir -p "$W/runs" "$W/logs"
# Refuse partial existing outputs rather than overwrite or silently skip a cell.
for P in D H; do for S in 0 1 2 3 4; do test ! -e "$W/runs/${P}_NORT_s$S"; done; done
printf 'set\tseed\tjob\n' > "$W/submitted.tsv"
for P in D H; do
 for S in 0 1 2 3 4; do
  job=$(sbatch --parsable -J "rtnull-${P}-s${S}" -o "$W/logs/%x.%j.out" rtnull/rtnull_hela.sbatch "$P" "$S")
  printf '%s\t%s\t%s\n' "$P" "$S" "$job" | tee -a "$W/submitted.tsv"
 done
done
