#!/bin/bash
#SBATCH --job-name=afull-predict
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=384G
#SBATCH --time=12:00:00
set -euo pipefail
cd /nfs/research/juan/singj/opendia/astral-full-20261123
mkdir -p library
scratch=$(mktemp -d "${TMPDIR:-/tmp}/afull-predict-${SLURM_JOB_ID}-XXXXXX")
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv > library/gpu.txt
sha256sum /nfs/research/juan/singj/opendia/proteobench-astral-entrapment/transfer/full_fasta/ProteoBenchFASTA_Entrapment_Human_with_contaminants_entrapment_pep.fasta /nfs/research/juan/singj/opendia/entlib-hela-20261118/dialibgen-cuda-0.13.1.sif > library/generation_inputs.sha256
/usr/bin/time -v -o library/generation.time.txt singularity exec --nv --cleanenv --no-home --home "$scratch:/tmp" --env TMPDIR=/tmp --bind /nfs/research/juan/singj/opendia --bind "$scratch:/tmp" /nfs/research/juan/singj/opendia/entlib-hela-20261118/dialibgen-cuda-0.13.1.sif DIALibGen \
 -in /nfs/research/juan/singj/opendia/proteobench-astral-entrapment/transfer/full_fasta/ProteoBenchFASTA_Entrapment_Human_with_contaminants_entrapment_pep.fasta \
 -out /nfs/research/juan/singj/opendia/astral-full-20261123/library/targets20.pqp \
 -threads 16 -machine:device cuda -machine:require_gpu \
 -generation:instrument Astral -generation:nce 25 \
 -generation:enzyme 'no cleavage' -generation:missed_cleavages 0 \
 -generation:precursor_mz 300 1800 -generation:precursor_charges 1 2 3 4 \
 -generation:fragment_mz 200 1800 -generation:peptide_length 7 30 \
 -generation:fixed_modifications 'Carbamidomethyl (C)' -generation:variable_modifications \
 -generation:n_terminal_methionine_excision false -generation:fragments 6 20 \
 -generation:max_fragment_charge 2 -generation:decoys none > library/generation.log 2>&1
sha256sum library/targets20.pqp > library/targets20.sha256
printf 'generation_complete\n' > library/generation.complete
