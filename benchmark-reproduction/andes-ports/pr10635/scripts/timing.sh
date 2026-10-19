#!/bin/bash
# Interleaved timing: PR build, astral_B1 and hfx_A2, plain DB, isotope_mimic off/on, 3 rounds, 4 threads, no Percolator.
set -euo pipefail
export LD_LIBRARY_PATH=/opt/openms-pr/lib OPENMS_DATA_PATH=/opt/openms-pr/share/OpenMS
for ds in astral_B1 hfx_A2; do for r in 1 2 3; do for arm in off on; do
  out=$S/runs/timing/$ds/$arm/r$r; rm -rf $out
  sets=(--set "modifications:variable=Oxidation (M)|Deamidated (N)|Deamidated (Q)" --set modifications:variable_max_per_peptide=2)
  [ $arm = on ] && sets+=(--set annotate:isotope_mimic=true)
  python3 -I $S/bench/run_search.py --bin /opt/openms-pr/bin --jobs $W/jobs.json --dataset $ds \
    --mzml $W/broad_benchmark/data/$ds.mzML --db $S/db/hye_td.fasta --out $out "${sets[@]}"
  echo "$ds $arm r$r $(cat $out/search.seconds) $(grep 'Maximum resident' $out/search.log | awk '{print $NF}')"
  rm -f $out/search.idXML $out/search.pin
done; done; done
