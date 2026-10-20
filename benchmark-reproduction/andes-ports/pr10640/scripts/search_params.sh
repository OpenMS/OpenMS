# ProSE parameters: frozen protocol (#10364) Lumos CID TMT arm, except: K carries TMT or Carbamyl (urea; variable, as in the deposited PD search), up to 2 variable mods; report top 1 (pooled Percolator)
PROSE_PARAMS="-Search:enzyme Trypsin/P -Search:decoys auto -Search:decoy_prefix DECOY_ \
 -Search:precursor:mass_tolerance_lower 20 -Search:precursor:mass_tolerance_upper 20 -Search:precursor:mass_tolerance_unit ppm \
 -Search:precursor:min_charge 2 -Search:precursor:max_charge 5 -Search:precursor:isotope_error_min -2 -Search:precursor:isotope_error_max 1 \
 -Search:fragment:mass_tolerance 0.5 -Search:fragment:mass_tolerance_unit Da -Search:fragment:min_ion_index 2 \
 -Search:modifications:fixed 'Carbamidomethyl (C)' 'TMT6plex (N-term)' -Search:modifications:variable 'TMT6plex (K)' 'Carbamyl (K)' 'Oxidation (M)' \
 -Search:modifications:variable_max_per_peptide 2 -Search:peptide:min_size 7 -Search:peptide:max_size 40 -Search:peptide:missed_cleavages 2 \
 -Search:peptide:deduplicate true -Search:report:top_hits 1 -Search:scoring:max_candidates_per_spectrum 50 \
 -Search:FDR:PSM 0 -Search:FDR:protein 0 -Search:rescore false -Search:calibration:enabled false"
