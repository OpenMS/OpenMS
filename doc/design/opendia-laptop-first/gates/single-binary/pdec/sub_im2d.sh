#!/bin/bash
set -e
W=/nfs/research/juan/singj/opendia/entlib-im2d-20261119
mkdir -p $W/runs
cat > $W/run.sbatch <<'EOF'
#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -t 4:00:00
# PLAN_IM2D (frozen 38273668e7): one OpenDIAWorkflow run on the r~7 predicted library. Usage: run.sbatch D|H SEED cur|new
C=/nfs/research/juan/singj/opendia; W=$C/entlib-im2d-20261119; V=$C/entlib-hela-20261118
R=$1; S=$2; A=$3; O=$W/runs/${R}_${A}_s$S
F=""; [ "$A" = new ] && F="-im_peaks 8"
[ "$R" = D ] && IN=$V/data/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
[ "$R" = H ] && IN=$V/data/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
mkdir -p $O; SC=${TMPDIR:-/tmp}/im2d-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $C/entlib-pdsif-20261118/opendia-pd9-07131138b8.sif OpenDIAWorkflow -in $IN -lib $C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp \
  -out_dir $O -threads 8 -seed $S $F > $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log; exit $rc
EOF
export PATH=/ebi/slurm/codon/bin:$PATH
for S in 0 1 2 3 4; do for A in cur new; do
  sbatch -J im_D_${A}_$S -o $W/runs/slurm-%j.out $W/run.sbatch D $S $A
done; done
date -u +%Y-%m-%dT%H:%M:%SZ
