#!/bin/bash
set -e
W=/nfs/research/juan/singj/opendia/entlib-im2d-20261119
export PATH=/ebi/slurm/codon/bin:$PATH
for S in 0 1 2 3 4; do for A in cur new; do
  sbatch -J im_H_${A}_$S -o $W/runs/slurm-%j.out $W/run.sbatch H $S $A
done; done
date -u +%Y-%m-%dT%H:%M:%SZ
