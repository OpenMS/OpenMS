#!/bin/bash
# PLAN_TRPPM: submit r~7 D (or H) seeds 0-4 x {cur,new}. Run on codon as pst_prd. Usage: sub_trppm.sh D|H
set -e
C=/nfs/research/juan/singj/opendia; W=$C/entlib-trppm-20261120
S=$C/entlib-pdsif-20261118/opendia-pd10-0942ffbdf3.sif
echo "fca04be9f0e10b7886b7c45a453e0689eea3b843fe2653917bdd1a8752546e8d  $S" | sha256sum -c
(cd $W && sha256sum -c SHA256SUMS)
mkdir -p $W/runs
cat > $W/run.sbatch <<'EOF'
#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -t 4:00:00
# PLAN_TRPPM: one OpenDIAWorkflow run on the r~7 predicted library. Usage: run.sbatch D|H SEED cur|new
C=/nfs/research/juan/singj/opendia; W=$C/entlib-trppm-20261120; V=$C/entlib-hela-20261118
R=$1; S=$2; A=$3; O=$W/runs/${R}_${A}_s$S
F=""; [ "$A" = new ] && F="-trace_tolerance_ppm 10"
[ "$R" = D ] && IN=$V/data/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
[ "$R" = H ] && IN=$V/data/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
mkdir -p $O; SC=${TMPDIR:-/tmp}/trppm-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp \
  --env LD_LIBRARY_PATH=$W/lib:/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/onnxruntime-linux-x64-1.23.2/lib:/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib \
  --bind $C --bind $SC:/tmp $C/entlib-pdsif-20261118/opendia-pd10-0942ffbdf3.sif $W/bin/OpenDIAWorkflow \
  -in $IN -lib $C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp -out_dir $O -threads 8 -seed $S $F > $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log; exit $rc
EOF
export PATH=/ebi/slurm/codon/bin:$PATH
R=${1:-D}
for A in cur new; do for SEED in 0 1 2 3 4; do
  sbatch -J trp_${R}_${A}_$SEED -o $W/runs/slurm-%j.out $W/run.sbatch $R $SEED $A
done; done
