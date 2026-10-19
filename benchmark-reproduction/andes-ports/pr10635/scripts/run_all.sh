#!/bin/bash
# usage: run_all.sh <build: pr|pr2|base> <arm: off|on|annotate|remove|nodeam> <db: plain|entrap> [datasets...]
# Env: S (work dir), W (bundle root with broad_benchmark/data), PERC (percolator exe)
set -euo pipefail
build=$1; arm=$2; dbk=$3; shift 3
DS=${@:-velos_125_R1 velos_5000_R2 velos_25000_R3 hfx_A2 hfx_B1 hfx_B3 astral_A2 astral_B1 astral_B3 lumos_lfq_5192 lumos_lfq_5194 lumos_lfq_5199}
export LD_LIBRARY_PATH=/opt/openms-$build/lib OPENMS_DATA_PATH=/opt/openms-$build/share/OpenMS
for ds in $DS; do
  case $ds in velos*|lumos*) db=hy;; *) db=hye;; esac
  [ "$dbk" = entrap ] && fasta=$S/db/${db}_entrap_td.fasta || fasta=$S/db/${db}_td.fasta
  out=$S/runs/$build/$dbk/$arm/$ds
  [ -f $out/s137/target.tsv ] && continue
  rm -rf $out
  sets=()
  if [ "$arm" != nodeam ]; then
    sets+=(--set "modifications:variable=Oxidation (M)|Deamidated (N)|Deamidated (Q)" --set modifications:variable_max_per_peptide=2)
  fi
  # on: the first PR head (annotate:isotope_mimic, flag + Percolator feature); annotate/remove: precursor:isotope_mimic
  [ "$arm" = on ] && sets+=(--set annotate:isotope_mimic=true)
  [ "$arm" = annotate ] && sets+=(--set precursor:isotope_mimic=annotate)
  [ "$arm" = remove ] && sets+=(--set precursor:isotope_mimic=remove)
  python3 -I $S/bench/run_search.py --bin /opt/openms-$build/bin --jobs $W/jobs.json --dataset $ds \
    --mzml $W/broad_benchmark/data/$ds.mzML --db $fasta --out $out --percolator $PERC "${sets[@]}"
  echo "done $build $dbk $arm $ds"
done
