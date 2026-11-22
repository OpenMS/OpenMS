#!/bin/bash
# PRECHECK_GAP2: 12 runs (D/H x r0/r7 x def/V7M6/V8M8, seed 0, -features_out, sibling default), then part A
# (precheck_ceil0.py on the def tables) and part B (precheck_gap2_b.py once everything finished).
set -e
C=/nfs/research/juan/singj/opendia; W=$C/gap2-20261122
mkdir -p $W/logs; export PATH=/ebi/slurm/codon/bin:$PATH
cat > $W/run.sbatch <<'EOF'
#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=160G
#SBATCH -t 6:00:00
# Usage: run.sbatch D|H r0|r7 def|V7M6|V8M8
C=/nfs/research/juan/singj/opendia; W=$C/gap2-20261122; LS=$C/libsize-20261121; V=$C/entlib-hela-20261118/data
SIF=$C/opendia-pd10-sibdef.sif; R=$1; L=$2; A=$3
[ "$R" = D ] && IN=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
[ "$R" = H ] && IN=$V/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
[ "$L" = r0 ] && LIB=$LS/r0/lib.pqp
[ "$L" = r7 ] && LIB=$C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp
F=""; [ "$A" = V7M6 ] && F="-valley_fraction 0.7 -max_peaks 6"; [ "$A" = V8M8 ] && F="-valley_fraction 0.8 -max_peaks 8"
O=$W/${R}_${L}_$A; mkdir -p $O; SC=${TMPDIR:-/tmp}/g2-$SLURM_JOB_ID; mkdir -p $SC
/usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
  $SIF OpenDIAWorkflow -in $IN -lib $LIB -features_out -out_dir $O -threads 8 -seed 0 $F > $O.log 2>&1
rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log; exit $rc
EOF
cat > $W/ana.sbatch <<'EOF'
#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=160G
#SBATCH -t 8:00:00
# Usage: ana.sbatch D|H r0|r7
C=/nfs/research/juan/singj/opendia; W=$C/gap2-20261122; LS=$C/libsize-20261121; R=$1; L=$2
[ "$L" = r0 ] && LIB=$LS/r0/lib.pqp && DN=$LS/r0/diann/$R/report.parquet && RR=0
[ "$L" = r7 ] && LIB=$C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp && DN=$C/entlib-diann5-20261119/$R/report.parquet && RR=7.0539
export OMP_NUM_THREADS=16
$C/imfeat-20261120/venv/bin/python $W/precheck_ceil0.py $W/${R}_${L}_def/*.features.tsv $LIB $DN ${R}_$L $RR > $W/ana_${R}_$L.txt 2>&1
EOF
ALL=""
for R in D H; do for L in r0 r7; do
  J=$(sbatch --parsable -J g2_${R}_${L}_def -o $W/logs/%x-%j.out $W/run.sbatch $R $L def); ALL="$ALL:$J"
  A=$(sbatch --parsable -J g2_a_${R}_$L --dependency=afterok:$J -o $W/logs/%x-%j.out $W/ana.sbatch $R $L)
  for X in V7M6 V8M8; do J2=$(sbatch --parsable -J g2_${R}_${L}_$X -o $W/logs/%x-%j.out $W/run.sbatch $R $L $X); ALL="$ALL:$J2"; done
  echo "$R $L def $J ana $A"
done; done
echo "runs$ALL"
