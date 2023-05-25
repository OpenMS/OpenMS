inputs:
  include:
    doc: Inclusion list input file in FASTA or featureXML format.
    type: File?
  exclude:
    doc: Exclusion list input file in featureXML, idXML or FASTA format.
    type: File?
  out:
    doc: Output file (tab delimited csv file).
    type: string
  rt_model:
    doc: RTModel file used for the rt prediction of peptides in FASTA files.
    type: File?
  pt_model:
    doc: PTModel file used for the pt prediction of peptides in FASTA files (only needed for inclusion_strategy PreotinBased_LP).
    type: File?
  inclusion_charges:
    doc: List containing the charge states to be considered for the inclusion list compounds, space separated.
    type: long[]?
  inclusion_strategy:
    doc: strategy to be used for selection
    type: string?
  exclusion_charges:
    doc: List containing the charge states to be considered for the exclusion list compounds (for idXML and FASTA input), space separated.
    type: long[]?
  raw_data:
    doc: File containing the raw data (only needed for FeatureBased_LP).
    type: File?
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
  algorithm__InclusionExclusionList__missed_cleavages:
    doc: "Number of missed cleavages used for protein digestion.\n"
    type: long?
  algorithm__InclusionExclusionList__RT__unit:
    doc: Create lists with units as seconds instead of minutes
    type: string?
  algorithm__InclusionExclusionList__RT__use_relative:
    doc: Use relative RT window, which depends on RT of precursor.
    type: string?
  algorithm__InclusionExclusionList__RT__window_relative:
    doc: "[for RT:use_relative == true] The relative factor X for the RT exclusion window, e.g. the window is calculated as [rt - rt*X, rt + rt*X]."
    type: double?
  algorithm__InclusionExclusionList__RT__window_absolute:
    doc: "[for RT:use_relative == false] The absolute value X for the RT exclusion window in [sec], e.g. the window is calculated as [rt - X, rt + X]."
    type: double?
  algorithm__InclusionExclusionList__merge__mz_tol:
    doc: Two inclusion/exclusion windows are merged when they (almost) overlap in RT (see 'rt_tol') and are close in m/z by this tolerance. Unit of this is defined in 'mz_tol_unit'.
    type: double?
  algorithm__InclusionExclusionList__merge__mz_tol_unit:
    doc: Unit of 'mz_tol'
    type: string?
  algorithm__InclusionExclusionList__merge__rt_tol:
    doc: Maximal RT delta (in seconds) which would allow two windows in RT to overlap (which causes merging the windows). Two inclusion/exclusion windows are merged when they (almost) overlap in RT and are close in m/z by this tolerance (see 'mz_tol'). Unit of this param is [seconds].
    type: double?
  algorithm__PrecursorSelection__ms2_spectra_per_rt_bin:
    doc: Number of allowed MS/MS spectra in a retention time bin.
    type: long?
  algorithm__PrecursorSelection__exclude_overlapping_peaks:
    doc: If true, overlapping or nearby peaks (within 'min_mz_peak_distance') are excluded for selection.
    type: boolean?
  algorithm__PrecursorSelection__Exclusion__use_dynamic_exclusion:
    doc: If true dynamic exclusion is applied.
    type: boolean?
  algorithm__PrecursorSelection__Exclusion__exclusion_time:
    doc: The time (in seconds) a feature is excluded.
    type: double?
  algorithm__PrecursorSelection__ProteinBasedInclusion__max_list_size:
    doc: The maximal number of precursors in the inclusion list.
    type: long?
  algorithm__PrecursorSelection__ProteinBasedInclusion__rt__min_rt:
    doc: Minimal rt in seconds.
    type: double?
  algorithm__PrecursorSelection__ProteinBasedInclusion__rt__max_rt:
    doc: Maximal rt in seconds.
    type: double?
  algorithm__PrecursorSelection__ProteinBasedInclusion__rt__rt_step_size:
    doc: rt step size in seconds.
    type: double?
  algorithm__PrecursorSelection__ProteinBasedInclusion__rt__rt_window_size:
    doc: rt window size in seconds.
    type: long?
  algorithm__PrecursorSelection__ProteinBasedInclusion__thresholds__min_protein_id_probability:
    doc: Minimal protein probability for a protein to be considered identified.
    type: double?
  algorithm__PrecursorSelection__ProteinBasedInclusion__thresholds__min_pt_weight:
    doc: Minimal pt weight of a precursor
    type: double?
  algorithm__PrecursorSelection__ProteinBasedInclusion__thresholds__min_mz:
    doc: Minimal mz to be considered in protein based LP formulation.
    type: double?
  algorithm__PrecursorSelection__ProteinBasedInclusion__thresholds__max_mz:
    doc: Minimal mz to be considered in protein based LP formulation.
    type: double?
  algorithm__PrecursorSelection__ProteinBasedInclusion__thresholds__use_peptide_rule:
    doc: Use peptide rule instead of minimal protein id probability
    type: boolean?
  algorithm__PrecursorSelection__ProteinBasedInclusion__thresholds__min_peptide_ids:
    doc: If use_peptide_rule is true, this parameter sets the minimal number of peptide ids for a protein id
    type: long?
  algorithm__PrecursorSelection__ProteinBasedInclusion__thresholds__min_peptide_probability:
    doc: If use_peptide_rule is true, this parameter sets the minimal probability for a peptide to be safely identified
    type: double?
  algorithm__PrecursorSelection__feature_based__no_intensity_normalization:
    doc: Flag indicating if intensities shall be scaled to be in [0,1]. This is done for each feature separately, so that the feature's maximal intensity in a spectrum is set to 1.
    type: boolean?
  algorithm__PrecursorSelection__feature_based__max_number_precursors_per_feature:
    doc: The maximal number of precursors per feature.
    type: long?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: InclusionExclusionListCreator
doc: Creates inclusion and/or exclusion lists.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - InclusionExclusionListCreator
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json