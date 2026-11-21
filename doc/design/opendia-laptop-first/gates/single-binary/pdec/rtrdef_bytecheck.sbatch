#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /nfs/research/juan/singj/opendia/rtrefine-20261121/bytecheck-%j.out
# Default flip byte check (r~1.56 D s0): new default vs tested -rt_refine arm, -no_rt_refine vs previous default.
C=/nfs/research/juan/singj/opendia
W=$C/rtrefine-20261121
V=$C/entlib-hela-20261118
L=/homes/sing/hela_entrap_codon.predicted_decoys.pqp
S=$C/opendia-pd10-rtrdef.sif
A=${A:-new}
O=$W/bytecheck/D_${A}_s0
mkdir -p "$O"
IN=$V/data/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
F=""
[ "$A" = off ] && F="-no_rt_refine"
SC=${TMPDIR:-/tmp}/rtrbc-$SLURM_JOB_ID
mkdir -p "$SC"
singularity exec --cleanenv --no-home --home "$SC:/tmp" --env TMPDIR=/tmp --bind "$C" --bind /homes/sing --bind "$SC:/tmp" \
  "$S" OpenDIAWorkflow -in "$IN" -lib "$L" -out_dir "$O" -threads 8 -seed 0 -features_out $F > "$O.log" 2>&1
rc=$?
rm -rf "$SC"
echo "exit $rc" >> "$O.log"
exit $rc
