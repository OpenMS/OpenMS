inputs:
  in:
    doc: Input files separated by blank
    type: File[]
  tr:
    doc: transition file ('TraML','tsv','pqp')
    type: File
  tr_type:
    doc: "input file type -- default: determined from file extension or content\n"
    type: string?
  tr_irt:
    doc: transition file ('TraML')
    type: File?
  tr_irt_nonlinear:
    doc: additional nonlinear transition file ('TraML')
    type: File?
  rt_norm:
    doc: RT normalization file (how to map the RTs of this run to the ones stored in the library). If set, tr_irt may be omitted.
    type: File?
  swath_windows_file:
    doc: "Optional, tab-separated file containing the SWATH windows for extraction: lower_offset upper_offset. Note that the first line is a header and will be skipped."
    type: File?
  sort_swath_maps:
    doc: Sort input SWATH files when matching to SWATH windows from swath_windows_file
    type: boolean?
  enable_ms1:
    doc: Extract the precursor ion trace(s) and use for scoring if present
    type: string?
  enable_ipf:
    doc: Enable additional scoring of identification assays using IPF (see online documentation)
    type: string?
  out_features:
    doc: output file
    type: string
  out_tsv:
    doc: TSV output file (mProphet-compatible TSV file)
    type: string
  out_osw:
    doc: OSW output file (PyProphet-compatible SQLite file)
    type: string
  out_chrom:
    doc: Also output all computed chromatograms output in mzML (chrom.mzML) or sqMass (SQLite format)
    type: string
  out_qc:
    doc: Optional QC meta data (charge distribution in MS1). Only works with mzML input files.
    type: string
  min_upper_edge_dist:
    doc: Minimal distance to the upper edge of a Swath window to still consider a precursor, in Thomson
    type: double?
  sonar:
    doc: data is scanning SWATH data
    type: boolean?
  pasef:
    doc: data is PASEF data
    type: boolean?
  rt_extraction_window:
    doc: Only extract RT around this value (-1 means extract over the whole range, a value of 600 means to extract around +/- 300 s of the expected elution).
    type: double?
  extra_rt_extraction_window:
    doc: Output an XIC with a RT-window by this much larger (e.g. to visually inspect a larger area of the chromatogram)
    type: double?
  ion_mobility_window:
    doc: "Extraction window in ion mobility dimension (in 1/k0 or milliseconds depending on library). This is the full window size, e.g. a value of 10 milliseconds would extract 5 milliseconds on either side. -1 means extract over the whole range or ion mobility is not present. (Default for diaPASEF data: 0.06 1/k0)"
    type: double?
  mz_extraction_window:
    doc: Extraction window in Thomson or ppm (see mz_extraction_window_unit)
    type: double?
  mz_extraction_window_unit:
    doc: Unit for mz extraction
    type: string?
  mz_extraction_window_ms1:
    doc: Extraction window used in MS1 in Thomson or ppm (see mz_extraction_window_ms1_unit)
    type: double?
  mz_extraction_window_ms1_unit:
    doc: Unit of the MS1 m/z extraction window
    type: string?
  im_extraction_window_ms1:
    doc: Extraction window in ion mobility dimension for MS1 (in 1/k0 or milliseconds depending on library). -1 means this is not ion mobility data.
    type: double?
  use_ms1_ion_mobility:
    doc: Also perform precursor extraction using the same ion mobility window as for fragment ion extraction
    type: string?
  matching_window_only:
    doc: Assume the input data is targeted / PRM-like data with potentially overlapping DIA windows. Will only attempt to extract each assay from the *best* matching DIA window (instead of all matching windows).
    type: boolean?
  irt_mz_extraction_window:
    doc: Extraction window used for iRT and m/z correction in Thomson or ppm (see irt_mz_extraction_window_unit)
    type: double?
  irt_mz_extraction_window_unit:
    doc: Unit for mz extraction
    type: string?
  irt_im_extraction_window:
    doc: Ion mobility extraction window used for iRT (in 1/K0 or milliseconds depending on library). -1 means do not perform ion mobility calibration
    type: double?
  min_rsq:
    doc: Minimum r-squared of RT peptides regression
    type: double?
  min_coverage:
    doc: Minimum relative amount of RT peptides to keep
    type: double?
  split_file_input:
    doc: "The input files each contain one single SWATH (alternatively: all SWATH are in separate files)"
    type: boolean?
  use_elution_model_score:
    doc: Turn on elution model score (EMG fit to peak)
    type: boolean?
  readOptions:
    doc: Whether to run OpenSWATH directly on the input data, cache data to disk first or to perform a datareduction step first. If you choose cache, make sure to also set tempDirectory
    type: string?
  mz_correction_function:
    doc: Use the retention time normalization peptide MS2 masses to perform a mass correction (linear, weighted by intensity linear or quadratic) of all spectra.
    type: string?
  tempDirectory:
    doc: Temporary directory to store cached files for example
    type: string?
  extraction_function:
    doc: Function used to extract the signal
    type: string?
  batchSize:
    doc: The batch size of chromatograms to process (0 means to only have one batch, sensible values are around 250-1000)
    type: long?
  outer_loop_threads:
    doc: How many threads should be used for the outer loop (-1 use all threads, use 4 to analyze 4 SWATH windows in memory at once).
    type: long?
  ms1_isotopes:
    doc: The number of MS1 isotopes used for extraction
    type: long?
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
  Debugging__irt_mzml:
    doc: Chromatogram mzML containing the iRT peptides
    type: string
  Debugging__irt_trafo:
    doc: Transformation file for RT transform
    type: string
  Calibration__ms1_im_calibration:
    doc: Whether to use MS1 precursor data for the ion mobility calibration (default = false, uses MS2 / fragment ions for calibration)
    type: boolean?
  Calibration__im_correction_function:
    doc: Type of normalization function for IM calibration.
    type: string?
  Calibration__debug_im_file:
    doc: Debug file for Ion Mobility calibration.
    type: string?
  Calibration__debug_mz_file:
    doc: Debug file for m/z calibration.
    type: string?
  Library__retentionTimeInterpretation:
    doc: How to interpret the provided retention time (the retention time column can either be interpreted to be in iRT, minutes or seconds)
    type: string?
  Library__override_group_label_check:
    doc: Override an internal check that assures that all members of the same PeptideGroupLabel have the same PeptideSequence (this ensures that only different isotopic forms of the same peptide can be grouped together in the same label group). Only turn this off if you know what you are doing.
    type: boolean?
  Library__force_invalid_mods:
    doc: Force reading even if invalid modifications are encountered (OpenMS may not recognize the modification)
    type: boolean?
  RTNormalization__alignmentMethod:
    doc: "How to perform the alignment to the normalized RT space using anchor points. 'linear': perform linear regression (for few anchor points). 'interpolated': Interpolate between anchor points (for few, noise-free anchor points). 'lowess' Use local regression (for many, noisy anchor points). 'b_spline' use b splines for smoothing."
    type: string?
  RTNormalization__outlierMethod:
    doc: "Which outlier detection method to use (valid: 'iter_residual', 'iter_jackknife', 'ransac', 'none'). Iterative methods remove one outlier at a time. Jackknife approach optimizes for maximum r-squared improvement while 'iter_residual' removes the datapoint with the largest residual error (removal by residual is computationally cheaper, use this with lots of peptides)."
    type: string?
  RTNormalization__useIterativeChauvenet:
    doc: Whether to use Chauvenet's criterion when using iterative methods. This should be used if the algorithm removes too many datapoints but it may lead to true outliers being retained.
    type: boolean?
  RTNormalization__RANSACMaxIterations:
    doc: Maximum iterations for the RANSAC outlier detection algorithm.
    type: long?
  RTNormalization__RANSACMaxPercentRTThreshold:
    doc: Maximum threshold in RT dimension for the RANSAC outlier detection algorithm (in percent of the total gradient). Default is set to 3% which is around +/- 4 minutes on a 120 gradient.
    type: long?
  RTNormalization__RANSACSamplingSize:
    doc: Sampling size of data points per iteration for the RANSAC outlier detection algorithm.
    type: long?
  RTNormalization__estimateBestPeptides:
    doc: Whether the algorithms should try to choose the best peptides based on their peak shape for normalization. Use this option you do not expect all your peptides to be detected in a sample and too many 'bad' peptides enter the outlier removal step (e.g. due to them being endogenous peptides or using a less curated list of peptides).
    type: boolean?
  RTNormalization__InitialQualityCutoff:
    doc: The initial overall quality cutoff for a peak to be scored (range ca. -2 to 2)
    type: double?
  RTNormalization__OverallQualityCutoff:
    doc: The overall quality cutoff for a peak to go into the retention time estimation (range ca. 0 to 10)
    type: double?
  RTNormalization__NrRTBins:
    doc: Number of RT bins to use to compute coverage. This option should be used to ensure that there is a complete coverage of the RT space (this should detect cases where only a part of the RT gradient is actually covered by normalization peptides)
    type: long?
  RTNormalization__MinPeptidesPerBin:
    doc: Minimal number of peptides that are required for a bin to counted as 'covered'
    type: long?
  RTNormalization__MinBinsFilled:
    doc: Minimal number of bins required to be covered
    type: long?
  RTNormalization__lowess__span:
    doc: Span parameter for lowess
    type: double?
  RTNormalization__b_spline__num_nodes:
    doc: Number of nodes for b spline
    type: long?
  Scoring__stop_report_after_feature:
    doc: Stop reporting after feature (ordered by quality; -1 means do not stop).
    type: long?
  Scoring__rt_normalization_factor:
    doc: The normalized RT is expected to be between 0 and 1. If your normalized RT has a different range, pass this here (e.g. it goes from 0 to 100, set this value to 100)
    type: double?
  Scoring__quantification_cutoff:
    doc: Cutoff in m/z below which peaks should not be used for quantification any more
    type: double?
  Scoring__write_convex_hull:
    doc: Whether to write out all points of all features into the featureXML
    type: boolean?
  Scoring__spectrum_addition_method:
    doc: For spectrum addition, either use simple concatenation or use peak resampling
    type: string?
  Scoring__add_up_spectra:
    doc: Add up spectra around the peak apex (needs to be a non-even integer)
    type: long?
  Scoring__spacing_for_spectra_resampling:
    doc: If spectra are to be added, use this spacing to add them up
    type: double?
  Scoring__uis_threshold_sn:
    doc: S/N threshold to consider identification transition (set to -1 to consider all)
    type: long?
  Scoring__uis_threshold_peak_area:
    doc: Peak area threshold to consider identification transition (set to -1 to consider all)
    type: long?
  Scoring__scoring_model:
    doc: Scoring model to use
    type: string?
  Scoring__im_extra_drift:
    doc: Extra drift time to extract for IM scoring (as a fraction, e.g. 0.25 means 25% extra on each side)
    type: double?
  Scoring__strict:
    doc: Whether to error (true) or skip (false) if a transition in a transition group does not have a corresponding chromatogram.
    type: string?
  Scoring__TransitionGroupPicker__stop_after_feature:
    doc: Stop finding after feature (ordered by intensity; -1 means do not stop).
    type: long?
  Scoring__TransitionGroupPicker__min_peak_width:
    doc: Minimal peak width (s), discard all peaks below this value (-1 means no action).
    type: double?
  Scoring__TransitionGroupPicker__peak_integration:
    doc: Calculate the peak area and height either the smoothed or the raw chromatogram data.
    type: string?
  Scoring__TransitionGroupPicker__background_subtraction:
    doc: "Remove background from peak signal using estimated noise levels. The 'original' method is only provided for historical purposes, please use the 'exact' method and set parameters using the PeakIntegrator: settings. The same original or smoothed chromatogram specified by peak_integration will be used for background estimation."
    type: string?
  Scoring__TransitionGroupPicker__recalculate_peaks:
    doc: Tries to get better peak picking by looking at peak consistency of all picked peaks. Tries to use the consensus (median) peak border if the variation within the picked peaks is too large.
    type: string?
  Scoring__TransitionGroupPicker__use_precursors:
    doc: Use precursor chromatogram for peak picking (note that this may lead to precursor signal driving the peak picking)
    type: boolean?
  Scoring__TransitionGroupPicker__use_consensus:
    doc: Use consensus peak boundaries when computing transition group picking (if false, compute independent peak boundaries for each transition)
    type: string?
  Scoring__TransitionGroupPicker__recalculate_peaks_max_z:
    doc: Determines the maximal Z-Score (difference measured in standard deviations) that is considered too large for peak boundaries. If the Z-Score is above this value, the median is used for peak boundaries (default value 1.0).
    type: double?
  Scoring__TransitionGroupPicker__minimal_quality:
    doc: Only if compute_peak_quality is set, this parameter will not consider peaks below this quality threshold
    type: double?
  Scoring__TransitionGroupPicker__resample_boundary:
    doc: For computing peak quality, how many extra seconds should be sample left and right of the actual peak
    type: double?
  Scoring__TransitionGroupPicker__compute_peak_quality:
    doc: Tries to compute a quality value for each peakgroup and detect outlier transitions. The resulting score is centered around zero and values above 0 are generally good and below -1 or -2 are usually bad.
    type: boolean?
  Scoring__TransitionGroupPicker__compute_peak_shape_metrics:
    doc: Calculates various peak shape metrics (e.g., tailing) that can be used for downstream QC/QA.
    type: boolean?
  Scoring__TransitionGroupPicker__compute_total_mi:
    doc: Compute mutual information metrics for individual transitions that can be used for OpenSWATH/IPF scoring.
    type: boolean?
  Scoring__TransitionGroupPicker__boundary_selection_method:
    doc: Method to use when selecting the best boundaries for peaks.
    type: string?
  Scoring__TransitionGroupPicker__PeakPickerMRM__sgolay_frame_length:
    doc: "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added."
    type: long?
  Scoring__TransitionGroupPicker__PeakPickerMRM__sgolay_polynomial_order:
    doc: Order of the polynomial that is fitted.
    type: long?
  Scoring__TransitionGroupPicker__PeakPickerMRM__gauss_width:
    doc: Gaussian width in seconds, estimated peak size.
    type: double?
  Scoring__TransitionGroupPicker__PeakPickerMRM__use_gauss:
    doc: Use Gaussian filter for smoothing (alternative is Savitzky-Golay filter)
    type: string?
  Scoring__TransitionGroupPicker__PeakPickerMRM__peak_width:
    doc: Force a certain minimal peak_width on the data (e.g. extend the peak at least by this amount on both sides) in seconds. -1 turns this feature off.
    type: double?
  Scoring__TransitionGroupPicker__PeakPickerMRM__signal_to_noise:
    doc: Signal-to-noise threshold at which a peak will not be extended any more. Note that setting this too high (e.g. 1.0) can lead to peaks whose flanks are not fully captured.
    type: double?
  Scoring__TransitionGroupPicker__PeakPickerMRM__write_sn_log_messages:
    doc: Write out log messages of the signal-to-noise estimator in case of sparse windows or median in rightmost histogram bin
    type: boolean?
  Scoring__TransitionGroupPicker__PeakPickerMRM__remove_overlapping_peaks:
    doc: Try to remove overlapping peaks during peak picking
    type: string?
  Scoring__TransitionGroupPicker__PeakPickerMRM__method:
    doc: Which method to choose for chromatographic peak-picking (OpenSWATH legacy on raw data, corrected picking on smoothed chromatogram or Crawdad on smoothed chromatogram).
    type: string?
  Scoring__TransitionGroupPicker__PeakIntegrator__integration_type:
    doc: The integration technique to use in integratePeak() and estimateBackground() which uses either the summed intensity, integration by Simpson's rule or trapezoidal integration.
    type: string?
  Scoring__TransitionGroupPicker__PeakIntegrator__baseline_type:
    doc: The baseline type to use in estimateBackground() based on the peak boundaries. A rectangular baseline shape is computed based either on the minimal intensity of the peak boundaries, the maximum intensity or the average intensity (base_to_base).
    type: string?
  Scoring__TransitionGroupPicker__PeakIntegrator__fit_EMG:
    doc: Fit the chromatogram/spectrum to the EMG peak model.
    type: string?
  Scoring__DIAScoring__dia_extraction_window:
    doc: DIA extraction window in Th or ppm.
    type: double?
  Scoring__DIAScoring__dia_extraction_unit:
    doc: DIA extraction window unit
    type: string?
  Scoring__DIAScoring__dia_centroided:
    doc: Use centroided DIA data.
    type: boolean?
  Scoring__DIAScoring__dia_byseries_intensity_min:
    doc: DIA b/y series minimum intensity to consider.
    type: double?
  Scoring__DIAScoring__dia_byseries_ppm_diff:
    doc: DIA b/y series minimal difference in ppm to consider.
    type: double?
  Scoring__DIAScoring__dia_nr_isotopes:
    doc: DIA number of isotopes to consider.
    type: long?
  Scoring__DIAScoring__dia_nr_charges:
    doc: DIA number of charges to consider.
    type: long?
  Scoring__DIAScoring__peak_before_mono_max_ppm_diff:
    doc: DIA maximal difference in ppm to count a peak at lower m/z when searching for evidence that a peak might not be monoisotopic.
    type: double?
  Scoring__EMGScoring__max_iteration:
    doc: Maximum number of iterations using by Levenberg-Marquardt algorithm.
    type: long?
  Scoring__EMGScoring__init_mom:
    doc: Initialize parameters using method of moments estimators.
    type: boolean?
  Scoring__Scores__use_shape_score:
    doc: Use the shape score (this score measures the similarity in shape of the transitions using a cross-correlation)
    type: string?
  Scoring__Scores__use_coelution_score:
    doc: Use the coelution score (this score measures the similarity in coelution of the transitions using a cross-correlation)
    type: string?
  Scoring__Scores__use_rt_score:
    doc: Use the retention time score (this score measure the difference in retention time)
    type: string?
  Scoring__Scores__use_library_score:
    doc: Use the library score
    type: string?
  Scoring__Scores__use_intensity_score:
    doc: Use the intensity score
    type: string?
  Scoring__Scores__use_nr_peaks_score:
    doc: Use the number of peaks score
    type: string?
  Scoring__Scores__use_total_xic_score:
    doc: Use the total XIC score
    type: string?
  Scoring__Scores__use_total_mi_score:
    doc: Use the total MI score
    type: boolean?
  Scoring__Scores__use_sn_score:
    doc: Use the SN (signal to noise) score
    type: string?
  Scoring__Scores__use_mi_score:
    doc: Use the MI (mutual information) score
    type: string?
  Scoring__Scores__use_dia_scores:
    doc: Use the DIA (SWATH) scores. If turned off, will not use fragment ion spectra for scoring.
    type: string?
  Scoring__Scores__use_ms1_correlation:
    doc: Use the correlation scores with the MS1 elution profiles
    type: boolean?
  Scoring__Scores__use_sonar_scores:
    doc: Use the scores for SONAR scans (scanning swath)
    type: boolean?
  Scoring__Scores__use_ion_mobility_scores:
    doc: Use the scores for Ion Mobility scans
    type: boolean?
  Scoring__Scores__use_ms1_fullscan:
    doc: Use the full MS1 scan at the peak apex for scoring (ppm accuracy of precursor and isotopic pattern)
    type: boolean?
  Scoring__Scores__use_ms1_mi:
    doc: Use the MS1 MI score
    type: string?
  Scoring__Scores__use_uis_scores:
    doc: Use UIS scores for peptidoform identification
    type: boolean?
  Scoring__Scores__use_ionseries_scores:
    doc: Use MS2-level b/y ion-series scores for peptidoform identification
    type: string?
  Scoring__Scores__use_ms2_isotope_scores:
    doc: Use MS2-level isotope scores (pearson & manhattan) across product transitions (based on ID if annotated or averagine)
    type: string?
outputs:
  out_features:
    type: File
    outputBinding:
      glob: $(inputs.out_features)
  out_tsv:
    type: File
    outputBinding:
      glob: $(inputs.out_tsv)
  out_osw:
    type: File
    outputBinding:
      glob: $(inputs.out_osw)
  out_chrom:
    type: File
    outputBinding:
      glob: $(inputs.out_chrom)
  out_qc:
    type: File
    outputBinding:
      glob: $(inputs.out_qc)
  Debugging__irt_mzml:
    type: File
    outputBinding:
      glob: $(inputs.Debugging__irt_mzml)
  Debugging__irt_trafo:
    type: File
    outputBinding:
      glob: $(inputs.Debugging__irt_trafo)
label: OpenSwathWorkflow
doc: Complete workflow to run OpenSWATH
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - OpenSwathWorkflow
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json