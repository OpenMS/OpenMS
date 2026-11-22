#!/bin/bash
# PRECHECK_V7M6T: one whole-node (48 of 48 cores) Ice Lake job per cell, ABBA (def, V7M6, V7M6, def) after a cache warm-up.
set -e
C=/nfs/research/juan/singj/opendia; W=$C/v7m6t-20261122
mkdir -p $W/logs; export PATH=/ebi/slurm/codon/bin:$PATH
cat > $W/cell.sbatch <<'EOF'
#!/bin/bash
#SBATCH -c 48
#SBATCH -C icelake
#SBATCH --mem=400G
#SBATCH -t 12:00:00
# Usage: cell.sbatch D|H r0|r7
C=/nfs/research/juan/singj/opendia; W=$C/v7m6t-20261122; LS=$C/libsize-20261121; V=$C/entlib-hela-20261118/data
SIF=$C/opendia-pd10-best.sif; R=$1; L=$2
[ "$R" = D ] && IN=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
[ "$R" = H ] && IN=$V/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
[ "$L" = r0 ] && LIB=$LS/r0/lib.pqp
[ "$L" = r7 ] && LIB=$C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp
echo "node $(hostname) $(lscpu | grep 'Model name')"
find $IN -type f -exec cat {} + > /dev/null; cat $LIB $SIF > /dev/null
for step in def:1 V7M6:1 V7M6:2 def:2; do
  A=${step%:*}; P=${step#*:}; F=""; [ "$A" = V7M6 ] && F="-valley_fraction 0.7 -max_peaks 6"
  O=$W/${R}_${L}_${A}_$P; mkdir -p $O; SC=${TMPDIR:-/tmp}/v7-$SLURM_JOB_ID-$P$A; mkdir -p $SC
  /usr/bin/time -v singularity exec --cleanenv --no-home --home $SC:/tmp --env TMPDIR=/tmp --bind $C --bind $SC:/tmp \
    $SIF OpenDIAWorkflow -in $IN -lib $LIB -out_dir $O -threads 8 -seed 0 $F > $O.log 2>&1
  rc=$?; rm -rf $SC; echo "exit $rc" >> $O.log; echo "$step exit $rc $(date +%T)"
done
EOF
for R in D H; do for L in r0 r7; do sbatch -J v7t_${R}_$L -o $W/logs/%x-%j.out $W/cell.sbatch $R $L; done; done
