#!/bin/bash
# Interleaved timing: astral_B1 and hfx_A2, plain DB, option off vs ARM2, 3 rounds, 4 threads, no Percolator.
# BUILD=pr ARM2=on: first PR head (annotate:isotope_mimic); BUILD=pr2 ARM2=remove: precursor:isotope_mimic=remove
set -euo pipefail
BUILD=${BUILD:-pr}; ARM2=${ARM2:-on}
export LD_LIBRARY_PATH=/opt/openms-$BUILD/lib OPENMS_DATA_PATH=/opt/openms-$BUILD/share/OpenMS
for ds in astral_B1 hfx_A2; do for r in 1 2 3; do for arm in off $ARM2; do
  out=$S/runs/timing/$BUILD/$ds/$arm/r$r; rm -rf $out
  sets=(--set "modifications:variable=Oxidation (M)|Deamidated (N)|Deamidated (Q)" --set modifications:variable_max_per_peptide=2)
  [ $arm = on ] && sets+=(--set annotate:isotope_mimic=true)
  [ $arm = remove ] && sets+=(--set precursor:isotope_mimic=remove)
  python3 -I $S/bench/run_search.py --bin /opt/openms-$BUILD/bin --jobs $W/jobs.json --dataset $ds \
    --mzml $W/broad_benchmark/data/$ds.mzML --db $S/db/hye_td.fasta --out $out "${sets[@]}"
  echo "$BUILD $ds $arm r$r $(cat $out/search.seconds) $(grep 'Maximum resident' $out/search.log | awk '{print $NF}')"
  rm -f $out/search.idXML $out/search.pin
done; done; done
