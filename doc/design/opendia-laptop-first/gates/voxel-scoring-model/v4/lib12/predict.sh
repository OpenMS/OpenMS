#!/bin/bash
# 12-fragment prediction for the same 100K-peptide FASTA (DIALibGen 0.13.1 release build, generate mode only).
# Settings chosen on a 200-protein smoke test to match the original 6-fragment library (Astral, NCE 25: 94 % of the original
# fragments reproduced, RT identical). At most 14 threads (CPUs 0-13), as requested.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-lib12-20261111
echo predict > STATUS
/usr/bin/time -v taskset -c 0-13 ../tools/dialibgen-0.13.1/bin/DIALibGen -in /home/singjc/Documents/datasets/proteobench-astral-entrapment/100K_sampled_fasta/ProteoBenchFASTA_Entrapment_Human_with_contaminants_entrapment_pep_100k_randomly_sampled.fasta -out pred12.parquet -threads 14 -machine:device cpu \
  -generation:instrument Astral -generation:nce 25 -generation:precursor_mz 400 1000 -generation:fragment_mz 100 2000 \
  -generation:precursor_charges 1 2 3 4 5 -generation:variable_modifications 'Oxidation (M)' -generation:max_variable_modifications 1 \
  -generation:missed_cleavages 2 -generation:peptide_length 7 35 -generation:n_terminal_methionine_excision false \
  -generation:fragments 1 12 -generation:max_fragment_charge 2 > predict.log 2>&1
echo $? > EXIT; echo done > STATUS
