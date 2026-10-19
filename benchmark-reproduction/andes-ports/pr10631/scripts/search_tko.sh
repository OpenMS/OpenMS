#!/bin/bash
# ProSE search of the PXD017803 TKO runs (identical for both arms; the PR does not touch identification).
# usage: search_tko.sh <stem> <fragment_tol> <fragment_unit> <threads>
source $(dirname "$0")/env.sh pr
stem=$1; ftol=$2; funit=$3; thr=${4:-2}
D=${DATA:-/home/user/data}
ProSE -in $D/pxd017803/$stem.mzML -database $D/db/yeast_sp.fasta -out_idxml $D/search/$stem.idXML \
  -Search:enzyme Trypsin/P -Search:peptide:missed_cleavages 2 -Search:peptide:min_size 7 -Search:peptide:max_size 40 \
  -Search:precursor:min_charge 2 -Search:precursor:max_charge 5 \
  -Search:precursor:mass_tolerance_lower 20 -Search:precursor:mass_tolerance_upper 20 -Search:precursor:mass_tolerance_unit ppm \
  -Search:precursor:isotope_error_min -2 -Search:precursor:isotope_error_max 1 \
  -Search:fragment:mass_tolerance $ftol -Search:fragment:mass_tolerance_unit $funit \
  -Search:modifications:fixed 'Carbamidomethyl (C)' 'TMTpro (K)' 'TMTpro (N-term)' -Search:modifications:variable 'Oxidation (M)' \
  -Search:modifications:variable_max_per_peptide 1 -Search:scoring:max_candidates_per_spectrum 50 -Search:report:top_hits 10 \
  -Search:peptide:deduplicate true -Search:FDR:PSM 0.01 -threads $thr -summary_out $D/search/$stem.summary.yaml > $D/search/$stem.log 2>&1
echo "exit $?"
