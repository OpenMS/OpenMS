#!/bin/bash
# One ProSE search plus Percolator (seeds 1, 42, 137) for one arm and one file.
# usage: run_search.sh <openms prefix> <arm> <mzML> <fasta> <outdir> [extra ProSE -Search: args...]
set -euo pipefail
PREFIX=$1; ARM=$2; MZML=$3; DB=$4; OUT=$5; shift 5
PERC=${PERCOLATOR:-percolator}
mkdir -p "$OUT"
export LD_LIBRARY_PATH="$PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OPENMS_DATA_PATH="$PREFIX/share/OpenMS"
if [ "${SKIP_SEARCH:-0}" != 1 ]; then
/usr/bin/time -v -o "$OUT/time.txt" "$PREFIX/bin/ProSE" -in "$MZML" -database "$DB" \
  -out_idxml "$OUT/search.idXML" -out_pin "$OUT/search.pin.tsv" -threads 4 \
  -Search:enzyme Trypsin/P -Search:decoy_prefix DECOY_ \
  -Search:peptide:missed_cleavages 2 -Search:peptide:min_size 7 -Search:peptide:max_size 40 \
  -Search:precursor:min_charge 2 -Search:precursor:max_charge 5 \
  -Search:precursor:mass_tolerance_lower 20 -Search:precursor:mass_tolerance_upper 20 -Search:precursor:mass_tolerance_unit ppm \
  -Search:fragment:mass_tolerance 20 -Search:fragment:mass_tolerance_unit ppm \
  -Search:modifications:fixed "Carbamidomethyl (C)" -Search:modifications:variable "Oxidation (M)" \
  -Search:modifications:variable_max_per_peptide 1 \
  -Search:scoring:max_candidates_per_spectrum 50 -Search:report:top_hits 10 -Search:peptide:deduplicate true \
  -Search:calibration:enabled false -Search:rescore false -Search:FDR:PSM 0 -Search:FDR:protein 0 \
  "$@" > "$OUT/prose.log" 2>&1
fi
[ "${SKIP_PERC:-0}" = 1 ] && exit 0
for seed in 1 42 137; do
  mkdir -p "$OUT/s$seed"
  "$PERC" -Y -U --seed $seed --results-psms "$OUT/s$seed/target.tsv" --decoy-results-psms "$OUT/s$seed/decoy.tsv" \
    "$OUT/search.pin.tsv" > "$OUT/s$seed/percolator.log" 2>&1
done
