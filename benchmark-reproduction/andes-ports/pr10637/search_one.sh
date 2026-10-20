#!/bin/bash
# Identification of one run as in #10585: Comet -> PSMFeatureExtractor -> Percolator, quantms-style q <= 0.10 filter, PEP as main score.
# Usage: search_one.sh <bin_dir> <mzML> <td_fasta> <out_dir> <comet.exe> <percolator>
set -euo pipefail
B=$1; MZ=$2; DB=$3; OUT=$4; COMET=$5; PERC=$6
r=$(basename "$MZ" .mzML); mkdir -p "$OUT"; T=$OUT/tmp_$r; mkdir -p "$T"
$B/CometAdapter -in "$MZ" -out $T/$r.comet.idXML -database "$DB" -comet_executable "$COMET" \
  -precursor_mass_tolerance 10 -precursor_error_units ppm -isotope_error 0/1 \
  -fragment_mass_tolerance 0.02 -fragment_bin_offset 0.0 -instrument high_res \
  -enzyme Trypsin -missed_cleavages 2 -min_peptide_length 7 -max_peptide_length 40 \
  -fixed_modifications "Carbamidomethyl (C)" -variable_modifications "Oxidation (M)" "Acetyl (Protein N-term)" \
  -max_variable_mods_in_peptide 3 -precursor_charge 2:5 -num_hits 1 \
  -PeptideIndexing:decoy_string DECOY_ -PeptideIndexing:IL_equivalent -PeptideIndexing:unmatched_action warn -threads 4 -no_progress
$B/PSMFeatureExtractor -in $T/$r.comet.idXML -out $T/$r.feat.idXML -threads 4 -no_progress
rm -f $T/$r.comet.idXML
$B/PercolatorAdapter -in $T/$r.feat.idXML -out $T/$r.perc.idXML -percolator_executable "$PERC" \
  -score_type pep -post_processing_tdc -decoy_pattern DECOY_ -subset_max_train 300000 -seed 42 -threads 4 -no_progress
rm -f $T/$r.feat.idXML
$B/IDScoreSwitcher -in $T/$r.perc.idXML -out $T/$r.q.idXML -new_score MS:1001491 -new_score_orientation lower_better -new_score_type q-value -old_score PEP -no_progress
$B/IDFilter -in $T/$r.q.idXML -out $T/$r.qf.idXML -score:psm 0.10 -no_progress
$B/IDScoreSwitcher -in $T/$r.qf.idXML -out $OUT/$r.idXML -new_score PEP -new_score_orientation lower_better -new_score_type "Posterior Error Probability" -old_score q-value -no_progress
rm -rf $T
