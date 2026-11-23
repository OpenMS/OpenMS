#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 4:00:00
# PRECHECK_GBMFIX Astral arm. Usage: sbatch -J gfx-RR-ARM-sS -o logs/%x.%j.out gbmfix_astral.sbatch RR S ARM
C=/nfs/research/juan/singj/opendia; W=$C/gbmfix2-20261123
R=$1; S=$2; A=$3
case $A in
  G0) F="";;
  G1) F="-gbm_rounds 1";;
  GQ) F="-gbm_positive_q 0.005";;
  *) echo "unknown arm $A"; exit 2;;
esac
O=$W/runs/raw${R}_${A}_s$S; mkdir -p $O
SC=${TMPDIR:-/tmp}/gfx-$SLURM_JOB_ID; mkdir -p $SC
echo "node $(hostname) cpu $(lscpu | grep 'Model name' | sed 's/.*: *//') flags $F" > $O.log
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $C/opendia-pd10-gbmvar.sif \
  OpenDIA -in $C/proteobench-astral-entrapment/transfer/raw/LFQ_Astral_DIA_15min_50ng_Human_$R.raw \
  -lib $C/gap-20261114/home/lib/library.pqp -scoring_lib $C/rtwin2-20261116/lib/library12_full.pqp \
  -out_dir $O -threads 16 -seed $S $F >> $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log
if [ "$A" = G0 ]; then cmp $O/LFQ_Astral_DIA_15min_50ng_Human_$R.precursors.tsv $C/astral-gbmfix-20261123/runs/raw${R}_s$S/LFQ_Astral_DIA_15min_50ng_Human_$R.precursors.tsv && echo "GATE SAME" >> $O.log || echo "GATE DIFF" >> $O.log; fi
exit $rc
