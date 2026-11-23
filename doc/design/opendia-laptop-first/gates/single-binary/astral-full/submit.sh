#!/bin/bash
set -euo pipefail
export PATH=/ebi/slurm/codon/bin:$PATH
cd /nfs/research/juan/singj/opendia/astral-full-20261123
test ! -e job_ids.txt
sha256sum -c PRECHECK_ASTRAL_FULL.sha256
mkdir -p logs
sha256sum /nfs/research/juan/singj/opendia/opendia-pd10-fastint2.sif /nfs/research/juan/singj/opendia/entlib-hela-20261118/dialibgen-cuda-0.13.1.sif /nfs/research/juan/singj/opendia/proteobench-astral-entrapment/transfer/full_fasta/ProteoBenchFASTA_Entrapment_Human_with_contaminants_entrapment_pep.fasta /nfs/research/juan/singj/opendia/proteobench-astral-entrapment/transfer/raw/LFQ_Astral_DIA_15min_50ng_Human_0{1,2,3}.raw > execution_inputs.sha256
prediction=$(sbatch --parsable --output=logs/predict-%j.log predict_gpu.sbatch)
printf '%s\n' "$prediction" > job_ids.txt
preparation=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$prediction" --output=logs/prepare-%j.log prepare.sbatch)
printf '%s\n' "$preparation" >> job_ids.txt
search=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$preparation" --output=logs/search-%A_%a.log run.sbatch)
printf '%s\n' "$search" >> job_ids.txt
evaluation=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterany:$search" --output=logs/evaluate-%j.log evaluate.sbatch)
printf '%s\n' "$evaluation" >> job_ids.txt
printf 'prediction=%s preparation=%s search=%s evaluation=%s\n' "$prediction" "$preparation" "$search" "$evaluation"
