#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 3:00:00
# PRECHECK_IMFEAT step 4 (diagnostic): re-extract with the C++-feasible core histogram variant, then compare to band6.
set -e
W=/nfs/research/juan/singj/opendia/imfeat-20261120
export OMP_NUM_THREADS=16 FAINT2=$W/faint2 CAL=$W/ims06.calibration.json FEATS=$W/ims06.features.tsv NPROC=16
cd $W
mkdir -p run2 && cp -n imfeat.log run2/ 2>/dev/null || true
OMP_NUM_THREADS=1 $W/venv/bin/python precheck_imfeat_extract.py > imfeat_extract2.log 2>&1
$W/venv/bin/python precheck_imfeat_core.py > imfeat_core.log 2>&1
