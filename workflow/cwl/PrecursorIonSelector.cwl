inputs:
  in:
    doc: Input feature map file (featureXML)
    type: File
  out:
    doc: modified feature map
    type: string
  next_feat:
    doc: feature map (featureXML) file with the selected precursors
    type: string
  ids:
    doc: file containing results of identification
    type: File
  num_precursors:
    doc: number of precursors to be selected
    type: long?
  raw_data:
    doc: Input profile data.
    type: File?
  load_preprocessing:
    doc: The preprocessed db is loaded from file, not calculated.
    type: boolean?
  store_preprocessing:
    doc: The preprocessed db is stored.
    type: boolean?
  simulation:
    doc: Simulate the whole LC-MS/MS run.
    type: boolean?
  sim_results:
    doc: File containing the results of the simulation run
    type: string
  db_path:
    doc: db file
    type: File?
  rt_model:
    doc: SVM Model for RTPredict
    type: File?
  dt_model:
    doc: SVM Model for PTPredict
    type: File?
  solver:
    doc: LP solver type
    type: string?
  fixed_modifications:
    doc: the modifications i.e. Carboxymethyl (C)
    type: string[]?
  log:
    doc: Name of log file (created only when specified)
    type: string?
  debug:
    doc: Sets the debug level
    type: long?
  threads:
    doc: Sets the number of threads allowed to be used by the TOPP tool
    type: long?
  no_progress:
    doc: Disables progress logging to command line
    type: boolean?
  force:
    doc: Overrides tool-specific checks
    type: boolean?
  algorithm__type:
    doc: Strategy for precursor ion selection.
    type: string?
  algorithm__max_iteration:
    doc: Maximal number of iterations.
    type: long?
  algorithm__rt_bin_capacity:
    doc: Maximal number of precursors per rt bin.
    type: long?
  algorithm__step_size:
    doc: Maximal number of precursors per iteration.
    type: long?
  algorithm__peptide_min_prob:
    doc: Minimal peptide probability.
    type: double?
  algorithm__sequential_spectrum_order:
    doc: If true, precursors are selected sequentially with respect to their RT.
    type: boolean?
  algorithm__MIPFormulation__thresholds__min_protein_probability:
    doc: Minimal protein probability for a protein to be considered in the ILP
    type: double?
  algorithm__MIPFormulation__thresholds__min_protein_id_probability:
    doc: Minimal protein probability for a protein to be considered identified.
    type: double?
  algorithm__MIPFormulation__thresholds__min_pt_weight:
    doc: Minimal pt weight of a precursor
    type: double?
  algorithm__MIPFormulation__thresholds__min_mz:
    doc: Minimal mz to be considered in protein based LP formulation.
    type: double?
  algorithm__MIPFormulation__thresholds__max_mz:
    doc: Minimal mz to be considered in protein based LP formulation.
    type: double?
  algorithm__MIPFormulation__thresholds__min_pred_pep_prob:
    doc: Minimal predicted peptide probability of a precursor
    type: double?
  algorithm__MIPFormulation__thresholds__min_rt_weight:
    doc: Minimal rt weight of a precursor
    type: double?
  algorithm__MIPFormulation__thresholds__use_peptide_rule:
    doc: Use peptide rule instead of minimal protein id probability
    type: boolean?
  algorithm__MIPFormulation__thresholds__min_peptide_ids:
    doc: If use_peptide_rule is true, this parameter sets the minimal number of peptide ids for a protein id
    type: long?
  algorithm__MIPFormulation__thresholds__min_peptide_probability:
    doc: If use_peptide_rule is true, this parameter sets the minimal probability for a peptide to be safely identified
    type: double?
  algorithm__MIPFormulation__combined_ilp__k1:
    doc: "combined ilp: weight for z_i"
    type: double?
  algorithm__MIPFormulation__combined_ilp__k2:
    doc: "combined ilp: weight for x_j,s*int_j,s"
    type: double?
  algorithm__MIPFormulation__combined_ilp__k3:
    doc: "combined ilp: weight for -x_j,s*w_j,s"
    type: double?
  algorithm__MIPFormulation__combined_ilp__scale_matching_probs:
    doc: flag if detectability * rt_weight shall be scaled to cover all [0,1]
    type: string?
  algorithm__MIPFormulation__feature_based__no_intensity_normalization:
    doc: Flag indicating if intensities shall be scaled to be in [0,1]. This is done for each feature separately, so that the feature's maximal intensity in a spectrum is set to 1.
    type: boolean?
  algorithm__MIPFormulation__feature_based__max_number_precursors_per_feature:
    doc: The maximal number of precursors per feature.
    type: long?
  algorithm__Preprocessing__precursor_mass_tolerance:
    doc: Precursor mass tolerance which is used to query the peptide database for peptides
    type: double?
  algorithm__Preprocessing__precursor_mass_tolerance_unit:
    doc: Precursor mass tolerance unit.
    type: string?
  algorithm__Preprocessing__preprocessed_db_path:
    doc: Path where the preprocessed database should be stored
    type: string?
  algorithm__Preprocessing__preprocessed_db_pred_rt_path:
    doc: Path where the predicted rts of the preprocessed database should be stored
    type: string?
  algorithm__Preprocessing__preprocessed_db_pred_dt_path:
    doc: Path where the predicted rts of the preprocessed database should be stored
    type: string?
  algorithm__Preprocessing__max_peptides_per_run:
    doc: Number of peptides for that the pt and rt are parallelly predicted.
    type: long?
  algorithm__Preprocessing__missed_cleavages:
    doc: Number of allowed missed cleavages.
    type: long?
  algorithm__Preprocessing__taxonomy:
    doc: Taxonomy
    type: string?
  algorithm__Preprocessing__tmp_dir:
    doc: Absolute path to tmp data directory used to store files needed for rt and dt prediction.
    type: string?
  algorithm__Preprocessing__store_peptide_sequences:
    doc: Flag if peptide sequences should be stored.
    type: boolean?
  algorithm__Preprocessing__rt_settings__min_rt:
    doc: Minimal RT in the experiment (in seconds)
    type: double?
  algorithm__Preprocessing__rt_settings__max_rt:
    doc: Maximal RT in the experiment (in seconds)
    type: double?
  algorithm__Preprocessing__rt_settings__rt_step_size:
    doc: Time between two consecutive spectra (in seconds)
    type: double?
  algorithm__Preprocessing__rt_settings__gauss_mean:
    doc: mean of the gauss curve
    type: double?
  algorithm__Preprocessing__rt_settings__gauss_sigma:
    doc: std of the gauss curve
    type: double?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  next_feat:
    type: File
    outputBinding:
      glob: $(inputs.next_feat)
  sim_results:
    type: File
    outputBinding:
      glob: $(inputs.sim_results)
label: PrecursorIonSelector
doc: PrecursorIonSelector
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - PrecursorIonSelector
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json