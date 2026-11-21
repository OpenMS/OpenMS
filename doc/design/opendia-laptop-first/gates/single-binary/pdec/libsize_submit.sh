#!/bin/bash
# LIBSIZE (diagnostic): build nested r0/r1/r3 libraries, then run DIA-NN and OpenDIA (current default) on D and H in parallel.
set -e
C=/nfs/research/juan/singj/opendia
W=$C/libsize-20261121
PY=$C/imfeat-20261120/venv/bin/python
SIF=$C/opendia-pd10-rtrdef.sif
DNN=$C/andes-dia-20261116/in/diann-2.2.0.sif
V=$C/entlib-hela-20261118/data
echo "f334b5fc66c91bb311bc240e09074170f1e7686c95fe96e25315224ab56fa1b9  $SIF" | sha256sum -c
mkdir -p $W/logs
export PATH=/ebi/slurm/codon/bin:$PATH

cat > $W/build.sbatch <<EOF
#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 3:00:00
#SBATCH -o $W/logs/build-%j.out
$PY $W/libsize.py build $W > $W/logs/build.log 2>&1
EOF

cat > $W/opendia.sbatch <<'EOF'
#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 4:00:00
# Usage: opendia.sbatch LEVEL RUN SEED
C=/nfs/research/juan/singj/opendia; W=$C/libsize-20261121; V=$C/entlib-hela-20261118/data
SIF=$C/opendia-pd10-rtrdef.sif
L=$1; R=$2; S=$3
[ "$R" = D ] && IN=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
[ "$R" = H ] && IN=$V/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
O=$W/$L/opendia/${R}_s$S; mkdir -p $O
SC=${TMPDIR:-/tmp}/ls-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $SIF OpenDIAWorkflow -in $IN -lib $W/$L/lib.pqp -out_dir $O -threads 8 -seed $S > $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log; exit $rc
EOF

cat > $W/diann.sbatch <<'EOF'
#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 4:00:00
# Usage: diann.sbatch LEVEL RUN   (flags as RESULT_DIANN_R5: diann_r5.sbatch)
C=/nfs/research/juan/singj/opendia; W=$C/libsize-20261121; V=$C/entlib-hela-20261118/data
L=$1; R=$2
[ "$R" = D ] && IN=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
[ "$R" = H ] && IN=$V/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
O=$W/$L/diann/$R; mkdir -p $O/tmp
/usr/bin/time -v singularity exec --cleanenv --bind $C $C/andes-dia-20261116/in/diann-2.2.0.sif /usr/diann-2.2.0/diann-linux \
  --f $IN --lib $W/$L/targets.tsv --threads 16 --verbose 1 --out $O/report.parquet --qvalue 1 --gen-decoys --unimod4 \
  --min-pr-charge 1 --max-pr-charge 5 --temp $O/tmp > $O/diann.log 2>&1
rc=$?; echo "exit $rc" >> $O/diann.log; rm -rf $O/tmp; exit $rc
EOF

B=$(sbatch --parsable -J ls_build $W/build.sbatch)
echo "build $B"
for L in r0 r1 r3; do
  for R in D H; do
    sbatch --parsable -J ls_dnn_${L}_$R --dependency=afterok:$B -o $W/logs/dnn-${L}-$R-%j.out $W/diann.sbatch $L $R
    for S in 0 1; do
      sbatch --parsable -J ls_od_${L}_${R}_$S --dependency=afterok:$B -o $W/logs/od-${L}-$R-$S-%j.out $W/opendia.sbatch $L $R $S
    done
  done
done
