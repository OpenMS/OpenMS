#!/bin/bash
#SBATCH --gres=gpu:a100:1
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH -t 02:00:00
# PLAN_DECOYGEN step 2 on codon: DIALibGen 0.13.1 prediction of the decoy peptides on one GPU (MS2, RT on CUDA) with 32 host
# threads. Same options as lib/predict.sh (p6: NCE 23, 100-1800, 5-6 fragments; p12: NCE 25, 100-2000, 1-12 fragments).
# Runs inside the PAIRLIB image (glibc 2.39, as the DIALibGen build), with the CUDA 12.8 runtime and cuDNN 9 from the NVIDIA
# wheels (cudalib/lib) on LD_LIBRARY_PATH. GPU and CPU results agree to float rounding (DIALibGen docs, tf32 off).
# Usage: sbatch -o logs/predict.%j.out predict_gpu.sbatch
W=/nfs/research/juan/singj/opendia/decoygen-20261113
SIF=/nfs/research/juan/singj/opendia/pairlib-20261113/container/opendia-pairlib-a8be5bed3c.sif
DLG=$W/tools/dialibgen-0.13.1/bin/DIALibGen
cd $W/lib
nvidia-smi --query-gpu=name,memory.total --format=csv > gpu.txt
run() {  # name, extra options
  out=$1; shift
  singularity exec --nv --cleanenv --bind $W:$W --env LD_LIBRARY_PATH=$W/cudalib/lib,OPENBLAS_NUM_THREADS=1 $SIF \
    /usr/bin/time -v $DLG -in decoy_peptides.fasta -out decoys_$out.parquet -threads 32 -machine:device cuda -machine:require_gpu \
    -generation:instrument Astral -generation:enzyme "no cleavage" -generation:precursor_mz 300 1200 -generation:precursor_charges 1 2 3 4 5 \
    -generation:variable_modifications 'Oxidation (M)' -generation:max_variable_modifications 1 -generation:peptide_length 7 35 \
    -generation:n_terminal_methionine_excision false -generation:max_fragment_charge 2 "$@" > $out.log 2>&1
  echo $? > ${out^^}_EXIT
}
run p6 -generation:nce 23 -generation:fragment_mz 100 1800 -generation:fragments 5 6
run p12 -generation:nce 25 -generation:fragment_mz 100 2000 -generation:fragments 1 12
sha256sum decoy_map.tsv decoy_peptides.fasta decoys_p6.parquet decoys_p12.parquet > predict.sha256
