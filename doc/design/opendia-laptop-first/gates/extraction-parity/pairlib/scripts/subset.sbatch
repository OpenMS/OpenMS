#!/bin/bash
#SBATCH -J pairlib-subset
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH -c 2
#SBATCH -o subset.%j.out
set -euo pipefail
module load zstd/1.5.7
T=/nfs/research/juan/singj/opendia/proteobench-astral-entrapment/transfer
W=/nfs/research/juan/singj/opendia/pairlib-20261113
cd $W
sha256sum -c $T/full_fasta/ProteoBenchFASTA_Entrapment_Human_with_contaminants_entrapment.decoys.tsv.zst.sha256 2>&1 | tail -1 || true
zstd -dc $T/full_fasta/ProteoBenchFASTA_Entrapment_Human_with_contaminants_entrapment.decoys.tsv.zst | \
  python3 scripts/subset_pairlib.py $T/100K_sampled_fasta/ProteoBenchFASTA_Entrapment_Human_with_contaminants_entrapment_pep_100k_randomly_sampled.fasta \
  $T/full_fasta/ProteoBenchFASTA_Entrapment_Human_with_contaminants_entrapment_pep.fasta pairlib.decoys.tsv pairlib.fasta > subset_summary.txt
zstd -q -T2 -f pairlib.decoys.tsv -o pairlib.decoys.tsv.zst
sha256sum pairlib.decoys.tsv pairlib.decoys.tsv.zst pairlib.fasta > pairlib.sha256
echo done
