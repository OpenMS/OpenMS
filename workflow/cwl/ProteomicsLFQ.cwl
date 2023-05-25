inputs:
  in:
    doc: Input files
    type: File[]
  ids:
    doc: "Identifications filtered at PSM level (e.g., q-value < 0.01).And annotated with PEP as main score.\nWe suggest using:\n1. PSMFeatureExtractor to annotate percolator features.\n2. PercolatorAdapter tool (score_type = 'q-value', -post-processing-tdc)\n3. IDFilter (pep:score = 0.05)\nTo obtain well calibrated PEPs and an initial reduction of PSMs\nID files must be provided in same order as spectra files."
    type: File[]
  design:
    doc: design file
    type: File?
  fasta:
    doc: fasta file
    type: File?
  out:
    doc: output mzTab file
    type: string
  out_msstats:
    doc: output MSstats input file
    type: string
  out_triqler:
    doc: output Triqler input file
    type: string
  out_cxml:
    doc: output consensusXML file
    type: string
  proteinFDR:
    doc: Protein FDR threshold (0.05=5%).
    type: double?
  picked_proteinFDR:
    doc: Use a picked protein FDR?
    type: boolean?
  psmFDR:
    doc: FDR threshold for sub-protein level (e.g. 0.05=5%). Use -FDR_type to choose the level. Cutoff is applied at the highest level. If Bayesian inference was chosen, it is equivalent with a peptide FDR
    type: double?
  FDR_type:
    doc: Sub-protein FDR level. PSM, PSM+peptide (best PSM q-value).
    type: string?
  protein_inference:
    doc: "Infer proteins:\naggregation  = aggregates all peptide scores across a protein (using the best score) \nbayesian     = computes a posterior probability for every protein based on a Bayesian network.\n               Note: 'bayesian' only uses and reports the best PSM per peptide."
    type: string?
  protein_quantification:
    doc: "Quantify proteins based on:\nunique_peptides = use peptides mapping to single proteins or a group of indistinguishable proteins(according to the set of experimentally identified peptides).\nstrictly_unique_peptides = use peptides mapping to a unique single protein only.\nshared_peptides = use shared peptides only for its best group (by inference score)"
    type: string?
  quantification_method:
    doc: "feature_intensity: MS1 signal.\nspectral_counting: PSM counts."
    type: string?
  targeted_only:
    doc: "true: Only ID based quantification.\nfalse: include unidentified features so they can be linked to identified ones (=match between runs)."
    type: boolean?
  transfer_ids:
    doc: "Requantification using mean of aligned RTs of a peptide feature.\nOnly applies to peptides that were quantified in more than 50% of all runs (of a fraction)."
    type: string?
  mass_recalibration:
    doc: Mass recalibration.
    type: boolean?
  alignment_order:
    doc: If star, aligns all maps to the reference with most IDs,if treeguided, calculates a guiding tree first.
    type: string?
  keep_feature_top_psm_only:
    doc: If false, also keeps lower ranked PSMs that have the top-scoring sequence as a candidate per feature in the same file.
    type: string?
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
  Seeding__intThreshold:
    doc: Peak intensity threshold applied in seed detection.
    type: double?
  Seeding__charge:
    doc: Charge range considered for untargeted feature seeds.
    type: string?
  Seeding__traceRTTolerance:
    doc: Combines all spectra in the tolerance window to stabilize identification of isotope patterns. Controls sensitivity (low value) vs. specificity (high value) of feature seeds.
    type: double?
  Centroiding__signal_to_noise:
    doc: Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)
    type: double?
  Centroiding__spacing_difference_gap:
    doc: The extension of a peak is stopped if the spacing between two subsequent data points exceeds 'spacing_difference_gap * min_spacing'. 'min_spacing' is the smaller of the two spacings from the peak apex to its two neighboring points. '0' to disable the constraint. Not applicable to chromatograms.
    type: double?
  Centroiding__spacing_difference:
    doc: Maximum allowed difference between points during peak extension, in multiples of the minimal difference between the peak apex and its two neighboring points. If this difference is exceeded a missing point is assumed (see parameter 'missing'). A higher value implies a less stringent peak definition, since individual signals within the peak are allowed to be further apart. '0' to disable the constraint. Not applicable to chromatograms.
    type: double?
  Centroiding__missing:
    doc: Maximum number of missing points allowed when extending a peak to the left or to the right. A missing data point occurs if the spacing between two subsequent data points exceeds 'spacing_difference * min_spacing'. 'min_spacing' is the smaller of the two spacings from the peak apex to its two neighboring points. Not applicable to chromatograms.
    type: long?
  Centroiding__ms_levels:
    doc: List of MS levels for which the peak picking is applied. If empty, auto mode is enabled, all peaks which aren't picked yet will get picked. Other scans are copied to the output without changes.
    type: long[]?
  Centroiding__report_FWHM:
    doc: Add metadata for FWHM (as floatDataArray named 'FWHM' or 'FWHM_ppm', depending on param 'report_FWHM_unit') for each picked peak.
    type: boolean?
  Centroiding__report_FWHM_unit:
    doc: Unit of FWHM. Either absolute in the unit of input, e.g. 'm/z' for spectra, or relative as ppm (only sensible for spectra, not chromatograms).
    type: string?
  Centroiding__SignalToNoise__max_intensity:
    doc: maximal intensity considered for histogram construction. By default, it will be calculated automatically (see auto_mode). Only provide this parameter if you know what you are doing (and change 'auto_mode' to '-1')! All intensities EQUAL/ABOVE 'max_intensity' will be added to the LAST histogram bin. If you choose 'max_intensity' too small, the noise estimate might be too small as well.  If chosen too big, the bins become quite large (which you could counter by increasing 'bin_count', which increases runtime). In general, the Median-S/N estimator is more robust to a manual max_intensity than the MeanIterative-S/N.
    type: long?
  Centroiding__SignalToNoise__auto_max_stdev_factor:
    doc: "parameter for 'max_intensity' estimation (if 'auto_mode' == 0): mean + 'auto_max_stdev_factor' * stdev"
    type: double?
  Centroiding__SignalToNoise__auto_max_percentile:
    doc: "parameter for 'max_intensity' estimation (if 'auto_mode' == 1): auto_max_percentile th percentile"
    type: long?
  Centroiding__SignalToNoise__auto_mode:
    doc: "method to use to determine maximal intensity: -1 --> use 'max_intensity'; 0 --> 'auto_max_stdev_factor' method (default); 1 --> 'auto_max_percentile' method"
    type: long?
  Centroiding__SignalToNoise__win_len:
    doc: window length in Thomson
    type: double?
  Centroiding__SignalToNoise__bin_count:
    doc: number of bins for intensity values
    type: long?
  Centroiding__SignalToNoise__min_required_elements:
    doc: minimum number of elements required in a window (otherwise it is considered sparse)
    type: long?
  Centroiding__SignalToNoise__noise_for_empty_window:
    doc: noise value used for sparse windows
    type: double?
  Centroiding__SignalToNoise__write_log_messages:
    doc: Write out log messages in case of sparse windows or median in rightmost histogram bin
    type: string?
  PeptideQuantification__candidates_out:
    doc: Optional output file with feature candidates.
    type: string
  PeptideQuantification__debug:
    doc: Debug level for feature detection.
    type: long?
  PeptideQuantification__quantify_decoys:
    doc: Whether decoy peptides should be quantified (true) or skipped (false).
    type: boolean?
  PeptideQuantification__min_psm_cutoff:
    doc: Minimum score for the best PSM of a spectrum to be used as seed. Use 'none' for no cutoff.
    type: string?
  PeptideQuantification__extract__batch_size:
    doc: Nr of peptides used in each batch of chromatogram extraction. Smaller values decrease memory usage but increase runtime.
    type: long?
  PeptideQuantification__extract__mz_window:
    doc: "m/z window size for chromatogram extraction (unit: ppm if 1 or greater, else Da/Th)"
    type: double?
  PeptideQuantification__extract__n_isotopes:
    doc: Number of isotopes to include in each peptide assay.
    type: long?
  PeptideQuantification__extract__isotope_pmin:
    doc: Minimum probability for an isotope to be included in the assay for a peptide. If set, this parameter takes precedence over 'extract:n_isotopes'.
    type: double?
  PeptideQuantification__extract__rt_quantile:
    doc: Quantile of the RT deviations between aligned internal and external IDs to use for scaling the RT extraction window
    type: double?
  PeptideQuantification__extract__rt_window:
    doc: RT window size (in sec.) for chromatogram extraction. If set, this parameter takes precedence over 'extract:rt_quantile'.
    type: double?
  PeptideQuantification__detect__min_peak_width:
    doc: Minimum elution peak width. Absolute value in seconds if 1 or greater, else relative to 'peak_width'.
    type: double?
  PeptideQuantification__detect__signal_to_noise:
    doc: Signal-to-noise threshold for OpenSWATH feature detection
    type: double?
  PeptideQuantification__detect__mapping_tolerance:
    doc: RT tolerance (plus/minus) for mapping peptide IDs to features. Absolute value in seconds if 1 or greater, else relative to the RT span of the feature.
    type: double?
  PeptideQuantification__svm__samples:
    doc: Number of observations to use for training ('0' for all)
    type: long?
  PeptideQuantification__svm__no_selection:
    doc: By default, roughly the same number of positive and negative observations, with the same intensity distribution, are selected for training. This aims to reduce biases, but also reduces the amount of training data. Set this flag to skip this procedure and consider all available observations (subject to 'svm:samples').
    type: boolean?
  PeptideQuantification__svm__xval_out:
    doc: "Output file: SVM cross-validation (parameter optimization) results"
    type: string
  PeptideQuantification__svm__kernel:
    doc: SVM kernel
    type: string?
  PeptideQuantification__svm__xval:
    doc: Number of partitions for cross-validation (parameter optimization)
    type: long?
  PeptideQuantification__svm__log2_C:
    doc: Values to try for the SVM parameter 'C' during parameter optimization. A value 'x' is used as 'C = 2^x'.
    type: double[]?
  PeptideQuantification__svm__log2_gamma:
    doc: Values to try for the SVM parameter 'gamma' during parameter optimization (RBF kernel only). A value 'x' is used as 'gamma = 2^x'.
    type: double[]?
  PeptideQuantification__svm__log2_p:
    doc: Values to try for the SVM parameter 'epsilon' during parameter optimization (epsilon-SVR only). A value 'x' is used as 'epsilon = 2^x'.
    type: double[]?
  PeptideQuantification__svm__epsilon:
    doc: Stopping criterion
    type: double?
  PeptideQuantification__svm__cache_size:
    doc: Size of the kernel cache (in MB)
    type: double?
  PeptideQuantification__svm__no_shrinking:
    doc: Disable the shrinking heuristics
    type: boolean?
  PeptideQuantification__svm__predictors:
    doc: Names of OpenSWATH scores to use as predictors for the SVM (comma-separated list)
    type: string?
  PeptideQuantification__svm__min_prob:
    doc: Minimum probability of correctness, as predicted by the SVM, required to retain a feature candidate
    type: double?
  PeptideQuantification__model__type:
    doc: Type of elution model to fit to features
    type: string?
  PeptideQuantification__model__add_zeros:
    doc: Add zero-intensity points outside the feature range to constrain the model fit. This parameter sets the weight given to these points during model fitting; '0' to disable.
    type: double?
  PeptideQuantification__model__unweighted_fit:
    doc: Suppress weighting of mass traces according to theoretical intensities when fitting elution models
    type: boolean?
  PeptideQuantification__model__no_imputation:
    doc: If fitting the elution model fails for a feature, set its intensity to zero instead of imputing a value from the initial intensity estimate
    type: boolean?
  PeptideQuantification__model__each_trace:
    doc: Fit elution model to each individual mass trace
    type: boolean?
  PeptideQuantification__model__check__min_area:
    doc: Lower bound for the area under the curve of a valid elution model
    type: double?
  PeptideQuantification__model__check__boundaries:
    doc: Time points corresponding to this fraction of the elution model height have to be within the data region used for model fitting
    type: double?
  PeptideQuantification__model__check__width:
    doc: Upper limit for acceptable widths of elution models (Gaussian or EGH), expressed in terms of modified (median-based) z-scores. '0' to disable. Not applied to individual mass traces (parameter 'each_trace').
    type: double?
  PeptideQuantification__model__check__asymmetry:
    doc: Upper limit for acceptable asymmetry of elution models (EGH only), expressed in terms of modified (median-based) z-scores. '0' to disable. Not applied to individual mass traces (parameter 'each_trace').
    type: double?
  PeptideQuantification__EMGScoring__max_iteration:
    doc: Maximum number of iterations for EMG fitting.
    type: long?
  PeptideQuantification__EMGScoring__init_mom:
    doc: Alternative initial parameters for fitting through method of moments.
    type: boolean?
  Alignment__model_type:
    doc: Options to control the modeling of retention time transformations from data
    type: string?
  Alignment__model__type:
    doc: Type of model
    type: string?
  Alignment__model__linear__symmetric_regression:
    doc: Perform linear regression on 'y - x' vs. 'y + x', instead of on 'y' vs. 'x'.
    type: boolean?
  Alignment__model__linear__x_weight:
    doc: Weight x values
    type: string?
  Alignment__model__linear__y_weight:
    doc: Weight y values
    type: string?
  Alignment__model__linear__x_datum_min:
    doc: Minimum x value
    type: double?
  Alignment__model__linear__x_datum_max:
    doc: Maximum x value
    type: double?
  Alignment__model__linear__y_datum_min:
    doc: Minimum y value
    type: double?
  Alignment__model__linear__y_datum_max:
    doc: Maximum y value
    type: double?
  Alignment__model__b_spline__wavelength:
    doc: Determines the amount of smoothing by setting the number of nodes for the B-spline. The number is chosen so that the spline approximates a low-pass filter with this cutoff wavelength. The wavelength is given in the same units as the data; a higher value means more smoothing. '0' sets the number of nodes to twice the number of input points.
    type: double?
  Alignment__model__b_spline__num_nodes:
    doc: Number of nodes for B-spline fitting. Overrides 'wavelength' if set (to two or greater). A lower value means more smoothing.
    type: long?
  Alignment__model__b_spline__extrapolate:
    doc: "Method to use for extrapolation beyond the original data range. 'linear': Linear extrapolation using the slope of the B-spline at the corresponding endpoint. 'b_spline': Use the B-spline (as for interpolation). 'constant': Use the constant value of the B-spline at the corresponding endpoint. 'global_linear': Use a linear fit through the data (which will most probably introduce discontinuities at the ends of the data range)."
    type: string?
  Alignment__model__b_spline__boundary_condition:
    doc: "Boundary condition at B-spline endpoints: 0 (value zero), 1 (first derivative zero) or 2 (second derivative zero)"
    type: long?
  Alignment__model__lowess__span:
    doc: Fraction of datapoints (f) to use for each local regression (determines the amount of smoothing). Choosing this parameter in the range .2 to .8 usually results in a good fit.
    type: double?
  Alignment__model__lowess__num_iterations:
    doc: Number of robustifying iterations for lowess fitting.
    type: long?
  Alignment__model__lowess__delta:
    doc: Nonnegative parameter which may be used to save computations (recommended value is 0.01 of the range of the input, e.g. for data ranging from 1000 seconds to 2000 seconds, it could be set to 10). Setting a negative value will automatically do this.
    type: double?
  Alignment__model__lowess__interpolation_type:
    doc: "Method to use for interpolation between datapoints computed by lowess. 'linear': Linear interpolation. 'cspline': Use the cubic spline for interpolation. 'akima': Use an akima spline for interpolation"
    type: string?
  Alignment__model__lowess__extrapolation_type:
    doc: "Method to use for extrapolation outside the data range. 'two-point-linear': Uses a line through the first and last point to extrapolate. 'four-point-linear': Uses a line through the first and second point to extrapolate in front and and a line through the last and second-to-last point in the end. 'global-linear': Uses a linear regression to fit a line through all data points and use it for interpolation."
    type: string?
  Alignment__model__interpolated__interpolation_type:
    doc: Type of interpolation to apply.
    type: string?
  Alignment__model__interpolated__extrapolation_type:
    doc: "Type of extrapolation to apply: two-point-linear: use the first and last data point to build a single linear model, four-point-linear: build two linear models on both ends using the first two / last two points, global-linear: use all points to build a single linear model. Note that global-linear may not be continuous at the border."
    type: string?
  Alignment__align_algorithm__score_type:
    doc: Name of the score type to use for ranking and filtering (.oms input only). If left empty, a score type is picked automatically.
    type: string?
  Alignment__align_algorithm__score_cutoff:
    doc: Use only IDs above a score cut-off (parameter 'min_score') for alignment?
    type: boolean?
  Alignment__align_algorithm__min_score:
    doc: "If 'score_cutoff' is 'true': Minimum score for an ID to be considered.\nUnless you have very few runs or identifications, increase this value to focus on more informative peptides."
    type: double?
  Alignment__align_algorithm__min_run_occur:
    doc: "Minimum number of runs (incl. reference, if any) in which a peptide must occur to be used for the alignment.\nUnless you have very few runs or identifications, increase this value to focus on more informative peptides."
    type: long?
  Alignment__align_algorithm__max_rt_shift:
    doc: "Maximum realistic RT difference for a peptide (median per run vs. reference). Peptides with higher shifts (outliers) are not used to compute the alignment.\nIf 0, no limit (disable filter); if > 1, the final value in seconds; if <= 1, taken as a fraction of the range of the reference RT scale."
    type: double?
  Alignment__align_algorithm__use_unassigned_peptides:
    doc: Should unassigned peptide identifications be used when computing an alignment of feature or consensus maps? If 'false', only peptide IDs assigned to features will be used.
    type: boolean?
  Alignment__align_algorithm__use_feature_rt:
    doc: "When aligning feature or consensus maps, don't use the retention time of a peptide identification directly; instead, use the retention time of the centroid of the feature (apex of the elution profile) that the peptide was matched to. If different identifications are matched to one feature, only the peptide closest to the centroid in RT is used.\nPrecludes 'use_unassigned_peptides'."
    type: string?
  Alignment__align_algorithm__use_adducts:
    doc: If IDs contain adducts, treat differently adducted variants of the same molecule as different.
    type: string?
  Linking__use_identifications:
    doc: Never link features that are annotated with different peptides (only the best hit per peptide identification is taken into account).
    type: string?
  Linking__nr_partitions:
    doc: How many partitions in m/z space should be used for the algorithm (more partitions means faster runtime and more memory efficient execution).
    type: long?
  Linking__min_nr_diffs_per_bin:
    doc: "If IDs are used: How many differences from matching IDs should be used to calculate a linking tolerance for unIDed features in an RT region. RT regions will be extended until that number is reached."
    type: long?
  Linking__min_IDscore_forTolCalc:
    doc: "If IDs are used: What is the minimum score of an ID to assume a reliable match for tolerance calculation. Check your current score type!"
    type: double?
  Linking__noID_penalty:
    doc: "If IDs are used: For the normalized distances, how high should the penalty for missing IDs be? 0 = no bias, 1 = IDs inside the max tolerances always preferred (even if much further away)."
    type: double?
  Linking__ignore_charge:
    doc: "false [default]: pairing requires equal charge state (or at least one unknown charge '0'); true: Pairing irrespective of charge state"
    type: boolean?
  Linking__ignore_adduct:
    doc: "true [default]: pairing requires equal adducts (or at least one without adduct annotation); true: Pairing irrespective of adducts"
    type: string?
  Linking__distance_RT__exponent:
    doc: Normalized RT differences ([0-1], relative to 'max_difference') are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)
    type: double?
  Linking__distance_RT__weight:
    doc: Final RT distances are weighted by this factor
    type: double?
  Linking__distance_MZ__max_difference:
    doc: Never pair features with larger m/z distance (unit defined by 'unit')
    type: double?
  Linking__distance_MZ__unit:
    doc: Unit of the 'max_difference' parameter
    type: string?
  Linking__distance_MZ__exponent:
    doc: Normalized ([0-1], relative to 'max_difference') m/z differences are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)
    type: double?
  Linking__distance_MZ__weight:
    doc: Final m/z distances are weighted by this factor
    type: double?
  Linking__distance_intensity__exponent:
    doc: Differences in relative intensity ([0-1]) are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)
    type: double?
  Linking__distance_intensity__weight:
    doc: Final intensity distances are weighted by this factor
    type: double?
  Linking__distance_intensity__log_transform:
    doc: Log-transform intensities? If disabled, d = |int_f2 - int_f1| / int_max. If enabled, d = |log(int_f2 + 1) - log(int_f1 + 1)| / log(int_max + 1))
    type: string?
  ProteinQuantification__method:
    doc: "- top - quantify based on three most abundant peptides (number can be changed in 'top').\n- iBAQ (intensity based absolute quantification), calculate the sum of all peptide peak intensities divided by the number of theoretically observable tryptic peptides (https://rdcu.be/cND1J). Warning: only consensusXML or featureXML input is allowed!"
    type: string?
  ProteinQuantification__best_charge_and_fraction:
    doc: "Distinguish between fraction and charge states of a peptide. For peptides, abundances will be reported separately for each fraction and charge;\nfor proteins, abundances will be computed based only on the most prevalent charge observed of each peptide (over all fractions).\nBy default, abundances are summed over all charge states."
    type: boolean?
  ProteinQuantification__top__N:
    doc: Calculate protein abundance from this number of proteotypic peptides (most abundant first; '0' for all)
    type: long?
  ProteinQuantification__top__aggregate:
    doc: Aggregation method used to compute protein abundances from peptide abundances
    type: string?
  ProteinQuantification__top__include_all:
    doc: Include results for proteins with fewer proteotypic peptides than indicated by 'N' (no effect if 'N' is 0 or 1)
    type: string?
  ProteinQuantification__consensus__normalize:
    doc: Scale peptide abundances so that medians of all samples are equal
    type: boolean?
  ProteinQuantification__consensus__fix_peptides:
    doc: "Use the same peptides for protein quantification across all samples.\nWith 'N 0',all peptides that occur in every sample are considered.\nOtherwise ('N'), the N peptides that occur in the most samples (independently of each other) are selected,\nbreaking ties by total abundance (there is no guarantee that the best co-ocurring peptides are chosen!)."
    type: boolean?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_msstats:
    type: File
    outputBinding:
      glob: $(inputs.out_msstats)
  out_triqler:
    type: File
    outputBinding:
      glob: $(inputs.out_triqler)
  out_cxml:
    type: File
    outputBinding:
      glob: $(inputs.out_cxml)
  PeptideQuantification__candidates_out:
    type: File
    outputBinding:
      glob: $(inputs.PeptideQuantification__candidates_out)
  PeptideQuantification__svm__xval_out:
    type: File
    outputBinding:
      glob: $(inputs.PeptideQuantification__svm__xval_out)
label: ProteomicsLFQ
doc: A standard proteomics LFQ pipeline.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - ProteomicsLFQ
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json