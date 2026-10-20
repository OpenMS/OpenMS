#!/bin/bash
# Comet + PSMFeatureExtractor + Percolator, then quantms-style filter (q <= 0.10, PEP as main score); #10585 recipe
set -e; source $(dirname $0)/env.sh
r=$1; S=$WORK/search
$BIN/CometAdapter -in $DATA/$r.mzML -out $S/$r.comet.idXML -database $WORK/db/HYE_td.fasta -comet_executable $TP/Comet/comet.exe \
  -precursor_mass_tolerance 10 -precursor_error_units ppm -isotope_error 0/1 \
  -fragment_mass_tolerance 0.02 -fragment_bin_offset 0.0 -instrument high_res \
  -enzyme Trypsin -missed_cleavages 2 -min_peptide_length 7 -max_peptide_length 40 \
  -fixed_modifications "Carbamidomethyl (C)" -variable_modifications "Oxidation (M)" "Acetyl (Protein N-term)" \
  -max_variable_mods_in_peptide 3 -precursor_charge 2:5 -num_hits 1 \
  -PeptideIndexing:decoy_string DECOY_ -PeptideIndexing:IL_equivalent -PeptideIndexing:unmatched_action warn -threads 4 -no_progress
$BIN/PSMFeatureExtractor -in $S/$r.comet.idXML -out $S/$r.feat.idXML -threads 4 -no_progress
rm $S/$r.comet.idXML
$BIN/PercolatorAdapter -in $S/$r.feat.idXML -out $S/$r.perc.idXML -percolator_executable $TP/Percolator/percolator \
  -score_type pep -post_processing_tdc -decoy_pattern DECOY_ -subset_max_train 300000 -seed 42 -threads 4 -no_progress
rm $S/$r.feat.idXML
$BIN/IDScoreSwitcher -in $S/$r.perc.idXML -out $S/$r.q.idXML -new_score MS:1001491 -new_score_orientation lower_better -new_score_type q-value -old_score PEP -no_progress
$BIN/IDFilter -in $S/$r.q.idXML -out $S/$r.qf.idXML -score:psm 0.10 -no_progress
$BIN/IDScoreSwitcher -in $S/$r.qf.idXML -out $S/$r.idXML -new_score PEP -new_score_orientation lower_better -new_score_type "Posterior Error Probability" -old_score q-value -no_progress
rm $S/$r.q.idXML $S/$r.qf.idXML
