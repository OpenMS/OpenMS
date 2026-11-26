#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 6:00:00
# PRECHECK_EVSHARE Astral cell. Usage: sbatch -J ea-RR-ARM-sS -o logs/%x.%j.out evs_astral.sbatch RR S ARM   (ARM: A | G0)
C=/nfs/research/juan/singj/opendia; W=$C/evshare-20261126
R=$1; S=$2; A=$3
case $A in
  G0) F="";;
  A) F="-evidence_share";;
  *) echo "unknown arm $A"; exit 2;;
esac
O=$W/runs/raw${R}_${A}_s$S; test ! -e $O; mkdir -p $O
SC=${TMPDIR:-/tmp}/ea-$SLURM_JOB_ID; mkdir -p $SC
echo "node $(hostname) cpu $(lscpu | grep 'Model name' | sed 's/.*: *//') flags $F" > $O.log
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $W/opendia-evshare-20261126.sif \
  OpenDIA -in $C/proteobench-astral-entrapment/transfer/raw/LFQ_Astral_DIA_15min_50ng_Human_$R.raw \
  -lib $C/gap-20261114/home/lib/library.pqp -scoring_lib $C/rtwin2-20261116/lib/library12_full.pqp \
  -out_dir $O -threads 16 -seed $S $F >> $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log
REF=$C/astral-gbmfix-20261123/runs/raw${R}_s$S
if [ "$A" = G0 ] && [ -d $REF ]; then for f in $(cd $REF && ls); do cmp -s $O/$f $REF/$f && echo "GATE SAME $f" >> $O.log || echo "GATE DIFF $f" >> $O.log; done; fi
exit $rc
