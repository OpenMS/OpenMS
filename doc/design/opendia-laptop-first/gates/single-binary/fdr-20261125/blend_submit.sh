#!/bin/bash
# PRECHECK_BLEND submission (run on codon as pst_prd, in the new work dir, only after the frozen precheck is pushed).
# Step "gate": the 3 no-flag gate cells. Step "main": the 25 B5 cells. Step "freshg0": the 22 remaining G0 cells (only if a gate differs).
set -euo pipefail
export PATH=/ebi/slurm/codon/bin:$PATH
C=/nfs/research/juan/singj/opendia; W=$C/blend-20261125
cd $W
sha256sum -c inputs.sha256
mkdir -p logs
case $1 in
  gate)
    test ! -e gate_job_ids.txt
    { sbatch --parsable -J bl-01-G0-s0 -o logs/%x.%j.out blend_astral.sbatch 01 0 G0
      sbatch --parsable -J bh-D-G0-s0 -o logs/%x.%j.out blend_hela.sbatch D 0 G0
      sbatch --parsable -J bh-H-G0-s0 -o logs/%x.%j.out blend_hela.sbatch H 0 G0; } > gate_job_ids.txt;;
  main)
    test ! -e main_job_ids.txt
    { for R in 01 02 03; do for S in 0 1 2 3 4; do sbatch --parsable -J bl-$R-B5-s$S -o logs/%x.%j.out blend_astral.sbatch $R $S B5; done; done
      for P in D H; do for S in 0 1 2 3 4; do sbatch --parsable -J bh-$P-B5-s$S -o logs/%x.%j.out blend_hela.sbatch $P $S B5; done; done; } > main_job_ids.txt;;
  freshg0)
    test ! -e freshg0_job_ids.txt
    { for R in 01 02 03; do for S in 0 1 2 3 4; do [ "$R$S" = 010 ] || sbatch --parsable -J bl-$R-G0-s$S -o logs/%x.%j.out blend_astral.sbatch $R $S G0; done; done
      for P in D H; do for S in 1 2 3 4; do sbatch --parsable -J bh-$P-G0-s$S -o logs/%x.%j.out blend_hela.sbatch $P $S G0; done; done; } > freshg0_job_ids.txt;;
  *) echo "usage: submit.sh gate|main|freshg0"; exit 2;;
esac
