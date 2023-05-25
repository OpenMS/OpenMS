inputs:
  in:
    doc: "input file "
    type: File
  out:
    doc: "output file "
    type: string
  out_plot:
    doc: txt file (if gnuplot is available, a corresponding PDF will be created as well.)
    type: string
  split_charge:
    doc: The search engine scores are split by charge if this flag is set. Thus, for each charge state a new model will be computed.
    type: boolean?
  top_hits_only:
    doc: If set only the top hits of every PeptideIdentification will be used
    type: boolean?
  fdr_for_targets_smaller:
    doc: Only used, when top_hits_only set. Additionally, target/decoy information should be available. The score_type must be q-value from an previous False Discovery Rate run.
    type: double?
  ignore_bad_data:
    doc: If set errors will be written but ignored. Useful for pipelines with many datasets where only a few are bad, but the pipeline should run through.
    type: boolean?
  prob_correct:
    doc: If set scores will be calculated as '1 - ErrorProbabilities' and can be interpreted as probabilities for correct identifications.
    type: boolean?
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
  fit_algorithm__number_of_bins:
    doc: Number of bins used for visualization. Only needed if each iteration step of the EM-Algorithm will be visualized
    type: long?
  fit_algorithm__incorrectly_assigned:
    doc: for 'Gumbel', the Gumbel distribution is used to plot incorrectly assigned sequences. For 'Gauss', the Gauss distribution is used.
    type: string?
  fit_algorithm__max_nr_iterations:
    doc: Bounds the number of iterations for the EM algorithm when convergence is slow.
    type: long?
  fit_algorithm__neg_log_delta:
    doc: The negative logarithm of the convergence threshold for the likelihood increase.
    type: long?
  fit_algorithm__outlier_handling:
    doc: "What to do with outliers:\n- ignore_iqr_outliers: ignore outliers outside of 3*IQR from Q1/Q3 for fitting\n- set_iqr_to_closest_valid: set IQR-based outliers to the last valid value for fitting\n- ignore_extreme_percentiles: ignore everything outside 99th and 1st percentile (also removes equal values like potential censored max values in XTandem)\n- none: do nothing"
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_plot:
    type: File
    outputBinding:
      glob: $(inputs.out_plot)
label: IDPosteriorErrorProbability
doc: Estimates probabilities for incorrectly assigned peptide sequences and a set of search engine scores using a mixture model.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - IDPosteriorErrorProbability
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json