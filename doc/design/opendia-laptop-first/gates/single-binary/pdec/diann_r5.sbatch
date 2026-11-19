#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=160G
#SBATCH -t 12:00:00
#SBATCH -o /nfs/research/juan/singj/opendia/entlib-diann5-20261119/diann.%A_%a.out
# DIA-NN 2.2 on the r~7 entrapment target library, its own decoys (--gen-decoys), flags as in bench_hela.sh.
# Array 1 = D (2737.d), 2 = H (22474.d). Entrapment labels are used only in the evaluation.
C=/nfs/research/juan/singj/opendia; W=$C/entlib-diann5-20261119; V=$C/entlib-hela-20261118/data
[ "$SLURM_ARRAY_TASK_ID" = 1 ] && R=D && IN=$V/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
[ "$SLURM_ARRAY_TASK_ID" = 2 ] && R=H && IN=$V/20201207_tims03_Evo03_PS_SA_HeLa_200ng_EvoSep_prot_high_speed_21min_8cm_S1-C8_1_22474.d
O=$W/$R; mkdir -p $O/tmp
/usr/bin/time -v singularity exec --cleanenv --bind $C $C/andes-dia-20261116/in/diann-2.2.0.sif /usr/diann-2.2.0/diann-linux \
  --f $IN --lib $W/targets_r5.tsv --threads 16 --verbose 1 --out $O/report.parquet --qvalue 1 --gen-decoys --unimod4 \
  --min-pr-charge 1 --max-pr-charge 5 --temp $O/tmp > $O/diann.log 2>&1
rc=$?; echo "exit $rc" >> $O/diann.log; rm -rf $O/tmp; exit $rc
