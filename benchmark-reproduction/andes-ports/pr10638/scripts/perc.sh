#!/bin/bash
# Percolator (seeds 1, 42, 137) on the -out_pin of one arm. The PSMs of co-isolated precursors share the SpecId of
# their spectrum (they differ in ExpMass, which Percolator's target-decoy competition groups by), so input.pin gets
# unique SpecIds (<SpecId>#<row>) in every arm alike; nothing else changes.
# usage: perc.sh <arm dir>
set -euo pipefail
OUT=$1; PERC=${PERCOLATOR:-percolator}
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print; next} $1=="DefaultDirection"{print; next} {$1=$1"#"NR; print}' "$OUT/search.pin.tsv" > "$OUT/input.pin"
for seed in 1 42 137; do
  mkdir -p "$OUT/s$seed"
  "$PERC" -Y -U --seed $seed --results-psms "$OUT/s$seed/target.tsv" --decoy-results-psms "$OUT/s$seed/decoy.tsv" \
    "$OUT/input.pin" > "$OUT/s$seed/percolator.log" 2>&1
done
