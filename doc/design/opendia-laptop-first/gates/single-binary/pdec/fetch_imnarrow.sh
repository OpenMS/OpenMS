#!/bin/bash
# PLAN_IMNARROW evaluation fetch: per run copy the run precursors.tsv (gz) into ~/.jcode/scratch/dia2/gap7/imnarrow/<R>_<arm>_s<seed>/
# Usage: fetch_extra.sh D|H
R=$1
L=/home/singjc/.jcode/scratch/dia2/gap7/imnarrow
mkdir -p $L
for A in cur new; do for S in 0 1 2 3 4; do
  O=$L/${R}_${A}_s$S; mkdir -p $O
  ssh -o BatchMode=yes -o ControlMaster=no -o ControlPath=none codon "cd /nfs/research/juan/singj/opendia/entlib-imnarrow-20261120/runs && tail -1 ${R}_${A}_s$S.log >&2 && f=\$(ls ${R}_${A}_s$S/*.precursors.tsv | grep -v '/precursors.tsv') && gzip -c \$f" > $O/run.precursors.tsv.gz 2> $O/status.txt
  echo "${R}_${A}_s$S $(grep -v bashrc $O/status.txt | tail -1) $(stat -c %s $O/run.precursors.tsv.gz)"
done; done
