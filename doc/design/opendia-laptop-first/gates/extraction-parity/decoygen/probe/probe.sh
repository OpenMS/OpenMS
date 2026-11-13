#!/bin/bash
# Probe: which DIALibGen setting reproduces the ProteoBench library intensities (label-free; library only).
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/decoygen-20261113/probe
awk '/^>/{n++} n<=400' ../../pairlib-20261113/lib/pairlib.fasta > small.fasta
DLG=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/tools/dialibgen-0.13.1/bin/DIALibGen
run() {
  "$DLG" -in small.fasta -out "p_$1_$2.parquet" -generation:instrument "$1" -generation:nce "$2" -generation:fragment_mz 100 2000 \
    -generation:fragments 1 12 -threads 8 -machine:device cpu -generation:precursor_mz 400 1000 -generation:precursor_charges 1 2 3 4 5 \
    -generation:variable_modifications 'Oxidation (M)' -generation:max_variable_modifications 1 -generation:missed_cleavages 2 \
    -generation:peptide_length 7 35 -generation:n_terminal_methionine_excision false -generation:max_fragment_charge 2 > "log_$1_$2.txt" 2>&1
  echo "$1 $2 exit $?"
}
run Astral 25; run Astral 30; run QE 27; run QE 30; run Lumos 30
