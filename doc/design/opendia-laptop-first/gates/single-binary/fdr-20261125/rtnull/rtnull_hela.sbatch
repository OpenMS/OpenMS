#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=200G
#SBATCH -t 16:00:00
# Frozen PRECHECK_RTNULL: one existing D/H raw, one original seed, no RT refinement.
set -euo pipefail
C=/nfs/research/juan/singj/opendia
W=$C/rtnull-20261126
SIF=$C/evshare-20261126/opendia-evshare-20261126.sif
LIB=$C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp
V=$C/entlib-hela-20261118/data
P=$1
S=$2
case $P in
 D) IN=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d;;
 H) IN=$V/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d;;
 *) exit 2;;
esac
case $S in 0|1|2|3|4) ;; *) exit 2;; esac
O=$W/runs/${P}_NORT_s$S
test ! -e "$O"
mkdir -p "$O"
SC=${TMPDIR:-/tmp}/rtnull-$SLURM_JOB_ID
mkdir -p "$SC"
{
 echo "job=$SLURM_JOB_ID node=$(hostname) set=$P seed=$S flags=-no_rt_refine"
 echo 'container_sha256=c43320a6e0d350efa29d80db2cafe2470bd026edd5a4fbafe6eba6f5792f191a'
} > "$O.log"
set +e
/usr/bin/time -v singularity exec --cleanenv --no-home --home "$SC:/tmp" --env TMPDIR=/tmp --bind "$C" --bind "$SC:/tmp" \
 "$SIF" OpenDIA -in "$IN" -lib "$LIB" -out_dir "$O" -threads 8 -seed "$S" -no_rt_refine >> "$O.log" 2>&1
rc=$?
set -e
printf 'exit %s\n' "$rc" >> "$O.log"
printf '%s\n' "$rc" > "$O.exit"
# Leave per-job scratch for inspection. No reference inputs or outputs are edited.
exit "$rc"
