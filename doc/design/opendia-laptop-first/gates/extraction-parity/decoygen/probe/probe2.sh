#!/bin/bash
# Probe 2: does DIALibGen top-6 within 100-1800 m/z reproduce the library's 6-fragment selection?
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/decoygen-20261113/probe
DLG=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/tools/dialibgen-0.13.1/bin/DIALibGen
run() {
  "$DLG" -in small.fasta -out "q_$1_$2_$3.parquet" -generation:instrument Astral -generation:nce "$1" -generation:fragment_mz "$2" "$3" \
    -generation:fragments 1 6 -threads 8 -machine:device cpu -generation:precursor_mz 400 1000 -generation:precursor_charges 1 2 3 4 5 \
    -generation:variable_modifications 'Oxidation (M)' -generation:max_variable_modifications 1 -generation:missed_cleavages 2 \
    -generation:peptide_length 7 35 -generation:n_terminal_methionine_excision false -generation:max_fragment_charge 2 > "logq_$1_$2_$3.txt" 2>&1
  echo "$1 $2 $3 exit $?"
}
for n in ${NCES:-25}; do run $n 100 1800; done
