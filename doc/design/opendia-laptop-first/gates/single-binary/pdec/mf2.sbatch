#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 3:00:00
#SBATCH -o /nfs/research/juan/singj/opendia/imfeat-20261120/mf2-%j.out
# PRECHECK_MF2 (diagnostic): ingredient ablation of the raw-point matched filter against the new default (36 + mt_).
set -e
W=/nfs/research/juan/singj/opendia/imfeat-20261120
PY=$W/venv/bin/python
cd $W
export FAINT2=$W/faint2 CAL=$W/ims06.calibration.json NPROC=16 CUR_FEATS=x CUR_PREC=x OMP_NUM_THREADS=1
STEP=none STEP2=extract $PY precheck_mf2.py > mf2_extract.log 2>&1
export OMP_NUM_THREADS=16
R=/nfs/research/juan/singj/opendia/mtcheck-20261120/runs/D_mt_s0/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737
export NEW_FEATS=$R.features.tsv NEW_PREC=$R.precursors.tsv
STEP=none STEP2=eval $PY precheck_mf2.py > mf2_eval.log 2>&1
echo done >> mf2_eval.log
