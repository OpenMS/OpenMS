#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 3:00:00
#SBATCH -o /nfs/research/juan/singj/opendia/imfeat-20261120/mf-%j.out
# PRECHECK_MF (diagnostic): matched-filter raw-point features on the faint2 sample, then the frozen gate.
set -e
W=/nfs/research/juan/singj/opendia/imfeat-20261120
PY=$W/venv/bin/python
export OMP_NUM_THREADS=1 FAINT2=$W/faint2 CAL=$W/ims06.calibration.json NPROC=16
export CUR_FEATS=/nfs/research/juan/singj/opendia/faintcpp-20261120/runs/D_cur_s0/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.features.tsv
export CUR_PREC=/nfs/research/juan/singj/opendia/faintcpp-20261120/runs/D_cur_s0/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.precursors.tsv
cd $W
STEP=extract $PY precheck_mf.py > mf_extract.log 2>&1
export OMP_NUM_THREADS=16
STEP=eval $PY precheck_mf.py > mf_eval.log 2>&1
echo done >> mf_eval.log
