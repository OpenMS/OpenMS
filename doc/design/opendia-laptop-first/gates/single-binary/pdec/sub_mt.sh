#!/bin/bash
# PLAN_MT (frozen 20 Nov 2026): submit r~7 D seeds 0-4 x {cur,mt}.
set -e
C=/nfs/research/juan/singj/opendia
W=$C/entlib-mt-20261120
S=$C/opendia-pd10-34c36861f2.sif
echo "1cad85f5605b3c4a62ee7d860594aaaddcdcd1dc9099bd930ad5446f1ad16d06  $S" | sha256sum -c
mkdir -p "$W/runs"
cat > "$W/run.sbatch" <<'EOF'
#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -t 4:00:00
# PLAN_MT (frozen 20 Nov 2026): one r~7 run. Usage: run.sbatch D|H SEED cur|mt
C=/nfs/research/juan/singj/opendia
W=$C/entlib-mt-20261120
V=$C/entlib-hela-20261118
SIF=$C/opendia-pd10-34c36861f2.sif
LIB=$C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp
R=$1
SEED=$2
ARM=$3
O=$W/runs/${R}_${ARM}_s${SEED}
F=""
[ "$ARM" = mt ] && F="-matched_filter"
[ "$R" = D ] && IN=$V/data/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
[ "$R" = H ] && IN=$V/data/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
mkdir -p "$O"
SC=${TMPDIR:-/tmp}/mt-$SLURM_JOB_ID
mkdir -p "$SC"
/usr/bin/time -v singularity exec --cleanenv --no-home --home "$SC:/tmp" --env TMPDIR=/tmp --bind "$C" --bind "$SC:/tmp" \
  "$SIF" OpenDIAWorkflow -in "$IN" -lib "$LIB" -out_dir "$O" -threads 8 -seed "$SEED" $F > "$O.log" 2>&1
rc=$?
rm -rf "$SC"
echo "exit $rc" >> "$O.log"
exit $rc
EOF
export PATH=/ebi/slurm/codon/bin:$PATH
R=${1:-D}
for ARM in cur mt; do
  for SEED in 0 1 2 3 4; do
    sbatch -J mt_${R}_${ARM}_${SEED} -o "$W/runs/slurm-%j.out" "$W/run.sbatch" "$R" "$SEED" "$ARM"
  done
done
