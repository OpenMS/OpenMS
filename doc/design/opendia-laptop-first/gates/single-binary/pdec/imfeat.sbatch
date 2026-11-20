#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 6:00:00
# PRECHECK_IMFEAT (diagnostic): mobility hand features + 1-D / 2-D trace CNN ceilings on the faint2 sample.
# Threads are fixed to the allocation (16): extraction pool 16 workers, torch + sklearn 16 threads.
set -e
W=/nfs/research/juan/singj/opendia/imfeat-20261120
export UV_CACHE_DIR=$W/uvcache
if [ ! -x $W/venv/bin/python ]; then
  ~/.local/bin/uv venv --python 3.12 $W/venv
  ~/.local/bin/uv pip install --python $W/venv/bin/python numpy pandas pyarrow scikit-learn scipy pyopenms
  ~/.local/bin/uv pip install --python $W/venv/bin/python torch --index-url https://download.pytorch.org/whl/cpu
fi
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
export FAINT2=$W/faint2 CAL=$W/ims06.calibration.json FEATS=$W/ims06.features.tsv NPROC=16 NTHREADS=16
cd $W
if [ ! -s $FAINT2/patches2d.npy ]; then
  OMP_NUM_THREADS=1 $W/venv/bin/python precheck_imfeat_extract.py > imfeat_extract.log 2>&1
fi
$W/venv/bin/python precheck_imfeat.py > imfeat.log 2>&1
echo done >> imfeat.log
