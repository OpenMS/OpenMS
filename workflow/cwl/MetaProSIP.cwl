inputs:
  in_mzML:
    doc: Centroided MS1 data
    type: File
  in_fasta:
    doc: Protein sequence database
    type: File
  out_csv:
    doc: Column separated file with feature fitting result.
    type: string
  out_peptide_centric_csv:
    doc: Column separated file with peptide centric result.
    type: string
  in_featureXML:
    doc: Feature data annotated with identifications (IDMapper)
    type: File
  r_executable:
    doc: "Path to the R executable (default: 'R')"
    type: File?
  mz_tolerance_ppm:
    doc: Tolerance in ppm
    type: double?
  rt_tolerance_s:
    doc: Tolerance window around feature rt for XIC extraction
    type: double?
  intensity_threshold:
    doc: Intensity threshold to collect peaks in the MS1 spectrum.
    type: double?
  correlation_threshold:
    doc: Correlation threshold for reporting a RIA
    type: double?
  xic_threshold:
    doc: Minimum correlation to mono-isotopic peak for retaining a higher isotopic peak. If featureXML from reference file is used it should be disabled (set to -1) as no mono-isotopic peak is expected to be present.
    type: double?
  decomposition_threshold:
    doc: Minimum R-squared of decomposition that must be achieved for a peptide to be reported.
    type: double?
  weight_merge_window:
    doc: Decomposition coefficients within +- this rate window will be combined
    type: double?
  min_correlation_distance_to_averagine:
    doc: Minimum difference in correlation between incorporation pattern and averagine pattern. Positive values filter all RIAs passing the correlation threshold but that also show a better correlation to an averagine peptide. Disabled for values <= -1
    type: double?
  pattern_15N_TIC_threshold:
    doc: The most intense peaks of the theoretical pattern contributing to at least this TIC fraction are taken into account.
    type: double?
  pattern_13C_TIC_threshold:
    doc: The most intense peaks of the theoretical pattern contributing to at least this TIC fraction are taken into account.
    type: double?
  pattern_2H_TIC_threshold:
    doc: The most intense peaks of the theoretical pattern contributing to at least this TIC fraction are taken into account.
    type: double?
  pattern_18O_TIC_threshold:
    doc: The most intense peaks of the theoretical pattern contributing to at least this TIC fraction are taken into account.
    type: double?
  heatmap_bins:
    doc: Number of RIA bins for heat map generation.
    type: long?
  plot_extension:
    doc: Extension used for plots (png|svg|pdf).
    type: string?
  qc_output_directory:
    doc: Output directory for the quality report
    type: string?
  labeling_element:
    doc: Which element (single letter code) is labeled.
    type: string?
  use_unassigned_ids:
    doc: Include identifications not assigned to a feature in pattern detection.
    type: boolean?
  use_averagine_ids:
    doc: Use averagine peptides as model to perform pattern detection on unidentified peptides.
    type: boolean?
  report_natural_peptides:
    doc: Whether purely natural peptides are reported in the quality report.
    type: boolean?
  filter_monoisotopic:
    doc: Try to filter out mono-isotopic patterns to improve detection of low RIA patterns
    type: boolean?
  cluster:
    doc: Perform grouping
    type: boolean?
  observed_peak_fraction:
    doc: Fraction of observed/expected peaks.
    type: double?
  min_consecutive_isotopes:
    doc: Minimum number of consecutive isotopic intensities needed.
    type: long?
  score_plot_yaxis_min:
    doc: The minimum value of the score axis. Values smaller than zero usually only make sense if the observed peak fraction is set to 0.
    type: double?
  collect_method:
    doc: How RIAs are collected.
    type: string?
  lowRIA_correlation_threshold:
    doc: Correlation threshold for reporting low RIA patterns. Disable and take correlation_threshold value for negative values.
    type: double?
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
outputs:
  out_csv:
    type: File
    outputBinding:
      glob: $(inputs.out_csv)
  out_peptide_centric_csv:
    type: File
    outputBinding:
      glob: $(inputs.out_peptide_centric_csv)
label: MetaProSIP
doc: Performs proteinSIP on peptide features for elemental flux analysis.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MetaProSIP
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json