#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 2:00:00
#SBATCH -o /nfs/research/juan/singj/opendia/entlib-diann5-20261119/conv.%j.out
# Targets-only r~7 library (the same targets.pqp every PLAN_PDEC / PLAN_ODP run used) as an OpenSWATH TSV for DIA-NN.
C=/nfs/research/juan/singj/opendia; W=$C/entlib-diann5-20261119
SC=${TMPDIR:-/tmp}/conv-$SLURM_JOB_ID; mkdir -p $SC
singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $C/entlib-pdsif-20261118/opendia-pd3-bc4e4f8d9b.sif TargetedFileConverter \
  -in $C/entlib-r5-20261118/targets.pqp -out $W/targets_r5.tsv -out_type tsv > $W/conv.log 2>&1
rc=$?; sha256sum $C/entlib-r5-20261118/targets.pqp $W/targets_r5.tsv >> $W/conv.log; wc -l $W/targets_r5.tsv >> $W/conv.log
rm -rf $SC; exit $rc
