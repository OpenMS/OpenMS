#!/bin/bash
# PRECHECK_EVSHARE submission (run on codon as pst_prd in the new work dir, only after the frozen precheck is pushed).
# Step "gate": 3 no-flag cells with existing G0 references. Step "main": 25 A cells and 10 fresh G0 cells on the
# new sets N1, N2. Step "freshg0": the 12 remaining G0 cells on raw03, D, H (only if a gate cell differs).
set -euo pipefail
export PATH=/ebi/slurm/codon/bin:$PATH
W=/nfs/research/juan/singj/opendia/evshare-20261126
cd $W
sha256sum -c inputs.sha256
(cd data && sha1sum -c pride.sha1)
mkdir -p logs
case $1 in
  gate)
    test ! -e gate_job_ids.txt
    { sbatch --parsable -J ea-03-G0-s0 -o logs/%x.%j.out evs_astral.sbatch 03 0 G0
      sbatch --parsable -J eh-D-G0-s0 -o logs/%x.%j.out evs_hela.sbatch D 0 G0
      sbatch --parsable -J eh-H-G0-s0 -o logs/%x.%j.out evs_hela.sbatch H 0 G0; } > gate_job_ids.txt;;
  main)
    test ! -e main_job_ids.txt
    { for S in 0 1 2 3 4; do sbatch --parsable -J ea-03-A-s$S -o logs/%x.%j.out evs_astral.sbatch 03 $S A; done
      for P in D H N1 N2; do for S in 0 1 2 3 4; do sbatch --parsable -J eh-$P-A-s$S -o logs/%x.%j.out evs_hela.sbatch $P $S A; done; done
      for P in N1 N2; do for S in 0 1 2 3 4; do sbatch --parsable -J eh-$P-G0-s$S -o logs/%x.%j.out evs_hela.sbatch $P $S G0; done; done; } > main_job_ids.txt;;
  freshg0)
    test ! -e freshg0_job_ids.txt
    { for S in 1 2 3 4; do sbatch --parsable -J ea-03-G0-s$S -o logs/%x.%j.out evs_astral.sbatch 03 $S G0; done
      for P in D H; do for S in 1 2 3 4; do sbatch --parsable -J eh-$P-G0-s$S -o logs/%x.%j.out evs_hela.sbatch $P $S G0; done; done; } > freshg0_job_ids.txt;;
  *) echo "usage: evs_submit.sh gate|main|freshg0"; exit 2;;
esac
