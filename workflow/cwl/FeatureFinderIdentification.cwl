inputs:
  in:
    doc: "Input file: LC-MS raw data"
    type: File
  id:
    doc: "Input file: Peptide identifications derived directly from 'in'"
    type: File
  id_ext:
    doc: "Input file: 'External' peptide identifications (e.g. from aligned runs)"
    type: File?
  out:
    doc: "Output file: Features"
    type: string
  lib_out:
    doc: "Output file: Assay library"
    type: string
  chrom_out:
    doc: "Output file: Chromatograms"
    type: string
  candidates_out:
    doc: "Output file: Feature candidates (before filtering and model fitting)"
    type: string
  candidates_in:
    doc: "Input file: Feature candidates from a previous run. If set, only feature classification and elution model fitting are carried out, if enabled. Many parameters are ignored."
    type: File?
  debug:
    doc: Sets the debug level
    type: long?
  quantify_decoys:
    doc: Whether decoy peptides should be quantified (true) or skipped (false).
    type: boolean?
  min_psm_cutoff:
    doc: Minimum score for the best PSM of a spectrum to be used as seed. Use 'none' for no cutoff.
    type: string?
  log:
    doc: Name of log file (created only when specified)
    type: string?
  threads:
    doc: Sets the number of threads allowed to be used by the TOPP tool
    type: long?
  no_progress:
    doc: Disables progress logging to command line
    type: boolean?
  force:
    doc: Overrides tool-specific checks
    type: boolean?
  extract__batch_size:
    doc: Nr of peptides used in each batch of chromatogram extraction. Smaller values decrease memory usage but increase runtime.
    type: long?
  extract__mz_window:
    doc: "m/z window size for chromatogram extraction (unit: ppm if 1 or greater, else Da/Th)"
    type: double?
  extract__n_isotopes:
    doc: Number of isotopes to include in each peptide assay.
    type: long?
  extract__isotope_pmin:
    doc: Minimum probability for an isotope to be included in the assay for a peptide. If set, this parameter takes precedence over 'extract:n_isotopes'.
    type: double?
  extract__rt_quantile:
    doc: Quantile of the RT deviations between aligned internal and external IDs to use for scaling the RT extraction window
    type: double?
  extract__rt_window:
    doc: RT window size (in sec.) for chromatogram extraction. If set, this parameter takes precedence over 'extract:rt_quantile'.
    type: double?
  detect__peak_width:
    doc: Expected elution peak width in seconds, for smoothing (Gauss filter). Also determines the RT extration window, unless set explicitly via 'extract:rt_window'.
    type: double?
  detect__min_peak_width:
    doc: Minimum elution peak width. Absolute value in seconds if 1 or greater, else relative to 'peak_width'.
    type: double?
  detect__signal_to_noise:
    doc: Signal-to-noise threshold for OpenSWATH feature detection
    type: double?
  detect__mapping_tolerance:
    doc: RT tolerance (plus/minus) for mapping peptide IDs to features. Absolute value in seconds if 1 or greater, else relative to the RT span of the feature.
    type: double?
  svm__samples:
    doc: Number of observations to use for training ('0' for all)
    type: long?
  svm__no_selection:
    doc: By default, roughly the same number of positive and negative observations, with the same intensity distribution, are selected for training. This aims to reduce biases, but also reduces the amount of training data. Set this flag to skip this procedure and consider all available observations (subject to 'svm:samples').
    type: boolean?
  svm__xval_out:
    doc: "Output file: SVM cross-validation (parameter optimization) results"
    type: string
  svm__kernel:
    doc: SVM kernel
    type: string?
  svm__xval:
    doc: Number of partitions for cross-validation (parameter optimization)
    type: long?
  svm__log2_C:
    doc: Values to try for the SVM parameter 'C' during parameter optimization. A value 'x' is used as 'C = 2^x'.
    type: double[]?
  svm__log2_gamma:
    doc: Values to try for the SVM parameter 'gamma' during parameter optimization (RBF kernel only). A value 'x' is used as 'gamma = 2^x'.
    type: double[]?
  svm__log2_p:
    doc: Values to try for the SVM parameter 'epsilon' during parameter optimization (epsilon-SVR only). A value 'x' is used as 'epsilon = 2^x'.
    type: double[]?
  svm__epsilon:
    doc: Stopping criterion
    type: double?
  svm__cache_size:
    doc: Size of the kernel cache (in MB)
    type: double?
  svm__no_shrinking:
    doc: Disable the shrinking heuristics
    type: boolean?
  svm__predictors:
    doc: Names of OpenSWATH scores to use as predictors for the SVM (comma-separated list)
    type: string?
  svm__min_prob:
    doc: Minimum probability of correctness, as predicted by the SVM, required to retain a feature candidate
    type: double?
  model__type:
    doc: Type of elution model to fit to features
    type: string?
  model__add_zeros:
    doc: Add zero-intensity points outside the feature range to constrain the model fit. This parameter sets the weight given to these points during model fitting; '0' to disable.
    type: double?
  model__unweighted_fit:
    doc: Suppress weighting of mass traces according to theoretical intensities when fitting elution models
    type: boolean?
  model__no_imputation:
    doc: If fitting the elution model fails for a feature, set its intensity to zero instead of imputing a value from the initial intensity estimate
    type: boolean?
  model__each_trace:
    doc: Fit elution model to each individual mass trace
    type: boolean?
  model__check__min_area:
    doc: Lower bound for the area under the curve of a valid elution model
    type: double?
  model__check__boundaries:
    doc: Time points corresponding to this fraction of the elution model height have to be within the data region used for model fitting
    type: double?
  model__check__width:
    doc: Upper limit for acceptable widths of elution models (Gaussian or EGH), expressed in terms of modified (median-based) z-scores. '0' to disable. Not applied to individual mass traces (parameter 'each_trace').
    type: double?
  model__check__asymmetry:
    doc: Upper limit for acceptable asymmetry of elution models (EGH only), expressed in terms of modified (median-based) z-scores. '0' to disable. Not applied to individual mass traces (parameter 'each_trace').
    type: double?
  EMGScoring__max_iteration:
    doc: Maximum number of iterations for EMG fitting.
    type: long?
  EMGScoring__init_mom:
    doc: Alternative initial parameters for fitting through method of moments.
    type: boolean?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  lib_out:
    type: File
    outputBinding:
      glob: $(inputs.lib_out)
  chrom_out:
    type: File
    outputBinding:
      glob: $(inputs.chrom_out)
  candidates_out:
    type: File
    outputBinding:
      glob: $(inputs.candidates_out)
  svm__xval_out:
    type: File
    outputBinding:
      glob: $(inputs.svm__xval_out)
label: FeatureFinderIdentification
doc: Detects features in MS1 data based on peptide identifications.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FeatureFinderIdentification
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json