#!/bin/bash
# ProSE search of one archive TMT run (HYE ProteoBench FASTA, sha256 d9ac434d...).
# usage: search_arch.sh <stem> <fragment_tol> <fragment_unit> <TMT6plex|TMTpro> <threads>
source $(dirname "$0")/env.sh pr
stem=$1; ftol=$2; funit=$3; tag=$4; thr=${5:-2}
D=${DATA:-/home/user/data}
ProSE -in $D/arch/$stem.mzML -database $D/db/hye.fasta -out_idxml $D/search/$stem.idXML \
  -Search:enzyme Trypsin/P -Search:peptide:missed_cleavages 2 -Search:peptide:min_size 7 -Search:peptide:max_size 40 \
  -Search:precursor:min_charge 2 -Search:precursor:max_charge 5 \
  -Search:precursor:mass_tolerance_lower 20 -Search:precursor:mass_tolerance_upper 20 -Search:precursor:mass_tolerance_unit ppm \
  -Search:precursor:isotope_error_min -2 -Search:precursor:isotope_error_max 1 \
  -Search:fragment:mass_tolerance $ftol -Search:fragment:mass_tolerance_unit $funit \
  -Search:modifications:fixed 'Carbamidomethyl (C)' "$tag (K)" "$tag (N-term)" -Search:modifications:variable 'Oxidation (M)' \
  -Search:modifications:variable_max_per_peptide 1 -Search:scoring:max_candidates_per_spectrum 50 -Search:report:top_hits 10 \
  -Search:peptide:deduplicate true -Search:FDR:PSM 0.01 -threads $thr -summary_out $D/search/$stem.summary.yaml > $D/search/$stem.log 2>&1
echo "exit $?"
