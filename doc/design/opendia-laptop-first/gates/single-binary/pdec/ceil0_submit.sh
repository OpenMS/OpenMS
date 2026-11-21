#!/bin/bash
# PRECHECK_CEIL0: 4 -features_out runs (D/H x r0/r7, seed 0, current default), then the ceiling analysis per run.
set -e
C=/nfs/research/juan/singj/opendia; W=$C/ceil0-20261121; LS=$C/libsize-20261121
mkdir -p $W/logs; export PATH=/ebi/slurm/codon/bin:$PATH
cat > $W/feat.sbatch <<'EOF'
#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -t 4:00:00
# Usage: feat.sbatch D|H r0|r7
C=/nfs/research/juan/singj/opendia; W=$C/ceil0-20261121; LS=$C/libsize-20261121; V=$C/entlib-hela-20261118/data
SIF=$C/opendia-pd10-rtrdef.sif; R=$1; L=$2
[ "$R" = D ] && IN=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d && N=20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737
[ "$R" = H ] && IN=$V/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d && N=20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474
[ "$L" = r0 ] && LIB=$LS/r0/lib.pqp && REF=$LS/r0/opendia/${R}_s0
[ "$L" = r7 ] && LIB=$C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp && REF=$C/entlib-rtr-20261121/runs/${R}_rtr_s0
O=$W/${R}_$L; mkdir -p $O; SC=${TMPDIR:-/tmp}/c0-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $SIF OpenDIAWorkflow -in $IN -lib $LIB -features_out -out_dir $O -threads 8 -seed 0 > $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log
cmp $O/$N.precursors.tsv $REF/$N.precursors.tsv && echo "SAME as $REF" >> $O.log || echo "DIFFERS from $REF" >> $O.log
exit $rc
EOF
cat > $W/ana.sbatch <<'EOF'
#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=160G
#SBATCH -t 8:00:00
# Usage: ana.sbatch D|H r0|r7
C=/nfs/research/juan/singj/opendia; W=$C/ceil0-20261121; LS=$C/libsize-20261121; R=$1; L=$2
[ "$L" = r0 ] && LIB=$LS/r0/lib.pqp && DN=$LS/r0/diann/$R/report.parquet && RR=0
[ "$L" = r7 ] && LIB=$C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp && DN=$C/entlib-diann5-20261119/$R/report.parquet && RR=7.0539
export OMP_NUM_THREADS=16
$C/imfeat-20261120/venv/bin/python $W/precheck_ceil0.py $W/${R}_$L/*.features.tsv $LIB $DN ${R}_$L $RR > $W/ana_${R}_$L.txt 2>&1
EOF
for R in D H; do for L in r0 r7; do
  J=$(sbatch --parsable -J c0_f_${R}_$L -o $W/logs/f-${R}-$L-%j.out $W/feat.sbatch $R $L)
  A=$(sbatch --parsable -J c0_a_${R}_$L --dependency=afterok:$J -o $W/logs/a-${R}-$L-%j.out $W/ana.sbatch $R $L)
  echo "$R $L feat $J ana $A"
done; done
