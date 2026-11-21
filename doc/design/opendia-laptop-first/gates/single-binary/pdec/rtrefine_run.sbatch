#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH -o /nfs/research/juan/singj/opendia/rtrefine-20261121/slurm-%j.out
# RT refine workflow check (r~1.56 D, no plan): one binary, arms def | rtr, seeds 0/1.
C=/nfs/research/juan/singj/opendia
W=$C/rtrefine-20261121
V=$C/entlib-hela-20261118
L=/homes/sing/hela_entrap_codon.predicted_decoys.pqp
S=$C/${SIF:-opendia-pd10-rtrefine.sif}
A=${A:-def}
SEED=${SEED:-0}
O=$W/runs/D_${A}_s${SEED}
mkdir -p "$O"
IN=$V/data/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
F=""
[ "$A" = rtr ] && F="-rt_refine"
SC=${TMPDIR:-/tmp}/rtr-$SLURM_JOB_ID
mkdir -p "$SC"
/usr/bin/time -v singularity exec --cleanenv --no-home --home "$SC:/tmp" --env TMPDIR=/tmp --bind "$C" --bind /homes/sing --bind "$SC:/tmp" \
  "$S" OpenDIAWorkflow -in "$IN" -lib "$L" -out_dir "$O" -threads 8 -seed "$SEED" -features_out $F > "$O.log" 2>&1
rc=$?
rm -rf "$SC"
echo "exit $rc" >> "$O.log"
exit $rc
