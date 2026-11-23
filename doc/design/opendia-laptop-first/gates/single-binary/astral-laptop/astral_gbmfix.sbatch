#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 4:00:00
# Astral re-run with the GBM final-round fix (OpenMS 88e5da3023). Laptop profile: 16 threads, 32 GB. Diagnostic, no pass rule.
# Same inputs as astral_laptop.sbatch (arm-E search library, library12_full scoring library, self-calibrating).
# Usage: sbatch -J agf-RR-sS -o logs/%x.%j.out astral_gbmfix.sbatch RR S
C=/nfs/research/juan/singj/opendia; W=$C/astral-gbmfix-20261123
R=$1; S=$2; O=$W/runs/raw${R}_s$S
mkdir -p $O
SC=${TMPDIR:-/tmp}/ast-$SLURM_JOB_ID; mkdir -p $SC
echo "node $(hostname) cpu $(lscpu | grep 'Model name' | sed 's/.*: *//')" > $O.log
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $C/opendia-pd10-gbmfix.sif \
  OpenDIA -in $C/proteobench-astral-entrapment/transfer/raw/LFQ_Astral_DIA_15min_50ng_Human_$R.raw \
  -lib $C/gap-20261114/home/lib/library.pqp -scoring_lib $C/rtwin2-20261116/lib/library12_full.pqp \
  -out_dir $O -threads 16 -seed $S >> $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log; exit $rc
