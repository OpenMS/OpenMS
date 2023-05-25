inputs:
  in:
    doc: Input files separated by blank
    type: File[]
  tr:
    doc: transition file with the RT peptides ('TraML' or 'csv')
    type: File
  out:
    doc: output file
    type: string
  rt_norm:
    doc: RT normalization file (how to map the RTs of this run to the ones stored in the library)
    type: File?
  min_rsq:
    doc: Minimum r-squared of RT peptides regression
    type: double?
  min_coverage:
    doc: Minimum relative amount of RT peptides to keep
    type: double?
  estimateBestPeptides:
    doc: Whether the algorithms should try to choose the best peptides based on their peak shape for normalization. Use this option you do not expect all your peptides to be detected in a sample and too many 'bad' peptides enter the outlier removal step (e.g. due to them being endogenous peptides or using a less curated list of peptides).
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
  algorithm__stop_report_after_feature:
    doc: Stop reporting after feature (ordered by quality; -1 means do not stop).
    type: long?
  algorithm__rt_extraction_window:
    doc: Only extract RT around this value (-1 means extract over the whole range, a value of 500 means to extract around +/- 500 s of the expected elution). For this to work, the TraML input file needs to contain normalized RT values.
    type: double?
  algorithm__rt_normalization_factor:
    doc: The normalized RT is expected to be between 0 and 1. If your normalized RT has a different range, pass this here (e.g. it goes from 0 to 100, set this value to 100)
    type: double?
  algorithm__quantification_cutoff:
    doc: Cutoff in m/z below which peaks should not be used for quantification any more
    type: double?
  algorithm__write_convex_hull:
    doc: Whether to write out all points of all features into the featureXML
    type: boolean?
  algorithm__spectrum_addition_method:
    doc: For spectrum addition, either use simple concatenation or use peak resampling
    type: string?
  algorithm__add_up_spectra:
    doc: Add up spectra around the peak apex (needs to be a non-even integer)
    type: long?
  algorithm__spacing_for_spectra_resampling:
    doc: If spectra are to be added, use this spacing to add them up
    type: double?
  algorithm__uis_threshold_sn:
    doc: S/N threshold to consider identification transition (set to -1 to consider all)
    type: long?
  algorithm__uis_threshold_peak_area:
    doc: Peak area threshold to consider identification transition (set to -1 to consider all)
    type: long?
  algorithm__scoring_model:
    doc: Scoring model to use
    type: string?
  algorithm__im_extra_drift:
    doc: Extra drift time to extract for IM scoring (as a fraction, e.g. 0.25 means 25% extra on each side)
    type: double?
  algorithm__strict:
    doc: Whether to error (true) or skip (false) if a transition in a transition group does not have a corresponding chromatogram.
    type: string?
  algorithm__TransitionGroupPicker__stop_after_feature:
    doc: Stop finding after feature (ordered by intensity; -1 means do not stop).
    type: long?
  algorithm__TransitionGroupPicker__stop_after_intensity_ratio:
    doc: Stop after reaching intensity ratio
    type: double?
  algorithm__TransitionGroupPicker__min_peak_width:
    doc: Minimal peak width (s), discard all peaks below this value (-1 means no action).
    type: double?
  algorithm__TransitionGroupPicker__peak_integration:
    doc: Calculate the peak area and height either the smoothed or the raw chromatogram data.
    type: string?
  algorithm__TransitionGroupPicker__background_subtraction:
    doc: "Remove background from peak signal using estimated noise levels. The 'original' method is only provided for historical purposes, please use the 'exact' method and set parameters using the PeakIntegrator: settings. The same original or smoothed chromatogram specified by peak_integration will be used for background estimation."
    type: string?
  algorithm__TransitionGroupPicker__recalculate_peaks:
    doc: Tries to get better peak picking by looking at peak consistency of all picked peaks. Tries to use the consensus (median) peak border if the variation within the picked peaks is too large.
    type: boolean?
  algorithm__TransitionGroupPicker__use_precursors:
    doc: Use precursor chromatogram for peak picking (note that this may lead to precursor signal driving the peak picking)
    type: boolean?
  algorithm__TransitionGroupPicker__use_consensus:
    doc: Use consensus peak boundaries when computing transition group picking (if false, compute independent peak boundaries for each transition)
    type: string?
  algorithm__TransitionGroupPicker__recalculate_peaks_max_z:
    doc: Determines the maximal Z-Score (difference measured in standard deviations) that is considered too large for peak boundaries. If the Z-Score is above this value, the median is used for peak boundaries (default value 1.0).
    type: double?
  algorithm__TransitionGroupPicker__minimal_quality:
    doc: Only if compute_peak_quality is set, this parameter will not consider peaks below this quality threshold
    type: double?
  algorithm__TransitionGroupPicker__resample_boundary:
    doc: For computing peak quality, how many extra seconds should be sample left and right of the actual peak
    type: double?
  algorithm__TransitionGroupPicker__compute_peak_quality:
    doc: Tries to compute a quality value for each peakgroup and detect outlier transitions. The resulting score is centered around zero and values above 0 are generally good and below -1 or -2 are usually bad.
    type: boolean?
  algorithm__TransitionGroupPicker__compute_peak_shape_metrics:
    doc: Calculates various peak shape metrics (e.g., tailing) that can be used for downstream QC/QA.
    type: boolean?
  algorithm__TransitionGroupPicker__compute_total_mi:
    doc: Compute mutual information metrics for individual transitions that can be used for OpenSWATH/IPF scoring.
    type: boolean?
  algorithm__TransitionGroupPicker__boundary_selection_method:
    doc: Method to use when selecting the best boundaries for peaks.
    type: string?
  algorithm__TransitionGroupPicker__PeakPickerMRM__sgolay_frame_length:
    doc: "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added."
    type: long?
  algorithm__TransitionGroupPicker__PeakPickerMRM__sgolay_polynomial_order:
    doc: Order of the polynomial that is fitted.
    type: long?
  algorithm__TransitionGroupPicker__PeakPickerMRM__gauss_width:
    doc: Gaussian width in seconds, estimated peak size.
    type: double?
  algorithm__TransitionGroupPicker__PeakPickerMRM__use_gauss:
    doc: Use Gaussian filter for smoothing (alternative is Savitzky-Golay filter)
    type: string?
  algorithm__TransitionGroupPicker__PeakPickerMRM__peak_width:
    doc: Force a certain minimal peak_width on the data (e.g. extend the peak at least by this amount on both sides) in seconds. -1 turns this feature off.
    type: double?
  algorithm__TransitionGroupPicker__PeakPickerMRM__signal_to_noise:
    doc: Signal-to-noise threshold at which a peak will not be extended any more. Note that setting this too high (e.g. 1.0) can lead to peaks whose flanks are not fully captured.
    type: double?
  algorithm__TransitionGroupPicker__PeakPickerMRM__sn_win_len:
    doc: Signal to noise window length.
    type: double?
  algorithm__TransitionGroupPicker__PeakPickerMRM__sn_bin_count:
    doc: Signal to noise bin count.
    type: long?
  algorithm__TransitionGroupPicker__PeakPickerMRM__write_sn_log_messages:
    doc: Write out log messages of the signal-to-noise estimator in case of sparse windows or median in rightmost histogram bin
    type: boolean?
  algorithm__TransitionGroupPicker__PeakPickerMRM__remove_overlapping_peaks:
    doc: Try to remove overlapping peaks during peak picking
    type: string?
  algorithm__TransitionGroupPicker__PeakPickerMRM__method:
    doc: Which method to choose for chromatographic peak-picking (OpenSWATH legacy on raw data, corrected picking on smoothed chromatogram or Crawdad on smoothed chromatogram).
    type: string?
  algorithm__TransitionGroupPicker__PeakIntegrator__integration_type:
    doc: The integration technique to use in integratePeak() and estimateBackground() which uses either the summed intensity, integration by Simpson's rule or trapezoidal integration.
    type: string?
  algorithm__TransitionGroupPicker__PeakIntegrator__baseline_type:
    doc: The baseline type to use in estimateBackground() based on the peak boundaries. A rectangular baseline shape is computed based either on the minimal intensity of the peak boundaries, the maximum intensity or the average intensity (base_to_base).
    type: string?
  algorithm__TransitionGroupPicker__PeakIntegrator__fit_EMG:
    doc: Fit the chromatogram/spectrum to the EMG peak model.
    type: string?
  algorithm__DIAScoring__dia_extraction_window:
    doc: DIA extraction window in Th or ppm.
    type: double?
  algorithm__DIAScoring__dia_extraction_unit:
    doc: DIA extraction window unit
    type: string?
  algorithm__DIAScoring__dia_centroided:
    doc: Use centroided DIA data.
    type: boolean?
  algorithm__DIAScoring__dia_byseries_intensity_min:
    doc: DIA b/y series minimum intensity to consider.
    type: double?
  algorithm__DIAScoring__dia_byseries_ppm_diff:
    doc: DIA b/y series minimal difference in ppm to consider.
    type: double?
  algorithm__DIAScoring__dia_nr_isotopes:
    doc: DIA number of isotopes to consider.
    type: long?
  algorithm__DIAScoring__dia_nr_charges:
    doc: DIA number of charges to consider.
    type: long?
  algorithm__DIAScoring__peak_before_mono_max_ppm_diff:
    doc: DIA maximal difference in ppm to count a peak at lower m/z when searching for evidence that a peak might not be monoisotopic.
    type: double?
  algorithm__EMGScoring__interpolation_step:
    doc: Sampling rate for the interpolation of the model function.
    type: double?
  algorithm__EMGScoring__tolerance_stdev_bounding_box:
    doc: Bounding box has range [minimim of data, maximum of data] enlarged by tolerance_stdev_bounding_box times the standard deviation of the data.
    type: double?
  algorithm__EMGScoring__max_iteration:
    doc: Maximum number of iterations using by Levenberg-Marquardt algorithm.
    type: long?
  algorithm__EMGScoring__init_mom:
    doc: Initialize parameters using method of moments estimators.
    type: boolean?
  algorithm__EMGScoring__statistics__mean:
    doc: Centroid position of the model.
    type: double?
  algorithm__EMGScoring__statistics__variance:
    doc: Variance of the model.
    type: double?
  algorithm__Scores__use_shape_score:
    doc: Use the shape score (this score measures the similarity in shape of the transitions using a cross-correlation)
    type: string?
  algorithm__Scores__use_coelution_score:
    doc: Use the coelution score (this score measures the similarity in coelution of the transitions using a cross-correlation)
    type: string?
  algorithm__Scores__use_rt_score:
    doc: Use the retention time score (this score measure the difference in retention time)
    type: string?
  algorithm__Scores__use_library_score:
    doc: Use the library score
    type: string?
  algorithm__Scores__use_elution_model_score:
    doc: Use the elution model (EMG) score (this score fits a gaussian model to the peak and checks the fit)
    type: string?
  algorithm__Scores__use_intensity_score:
    doc: Use the intensity score
    type: string?
  algorithm__Scores__use_nr_peaks_score:
    doc: Use the number of peaks score
    type: string?
  algorithm__Scores__use_total_xic_score:
    doc: Use the total XIC score
    type: string?
  algorithm__Scores__use_total_mi_score:
    doc: Use the total MI score
    type: boolean?
  algorithm__Scores__use_sn_score:
    doc: Use the SN (signal to noise) score
    type: string?
  algorithm__Scores__use_mi_score:
    doc: Use the MI (mutual information) score
    type: boolean?
  algorithm__Scores__use_dia_scores:
    doc: Use the DIA (SWATH) scores. If turned off, will not use fragment ion spectra for scoring.
    type: string?
  algorithm__Scores__use_ms1_correlation:
    doc: Use the correlation scores with the MS1 elution profiles
    type: boolean?
  algorithm__Scores__use_sonar_scores:
    doc: Use the scores for SONAR scans (scanning swath)
    type: boolean?
  algorithm__Scores__use_ion_mobility_scores:
    doc: Use the scores for Ion Mobility scans
    type: boolean?
  algorithm__Scores__use_ms1_fullscan:
    doc: Use the full MS1 scan at the peak apex for scoring (ppm accuracy of precursor and isotopic pattern)
    type: boolean?
  algorithm__Scores__use_ms1_mi:
    doc: Use the MS1 MI score
    type: boolean?
  algorithm__Scores__use_uis_scores:
    doc: Use UIS scores for peptidoform identification
    type: boolean?
  algorithm__Scores__use_ionseries_scores:
    doc: Use MS2-level b/y ion-series scores for peptidoform identification
    type: string?
  algorithm__Scores__use_ms2_isotope_scores:
    doc: Use MS2-level isotope scores (pearson & manhattan) across product transitions (based on ID if annotated or averagine)
    type: string?
  peptideEstimation__InitialQualityCutoff:
    doc: The initial overall quality cutoff for a peak to be scored (range ca. -2 to 2)
    type: double?
  peptideEstimation__OverallQualityCutoff:
    doc: The overall quality cutoff for a peak to go into the retention time estimation (range ca. 0 to 10)
    type: double?
  peptideEstimation__NrRTBins:
    doc: Number of RT bins to use to compute coverage. This option should be used to ensure that there is a complete coverage of the RT space (this should detect cases where only a part of the RT gradient is actually covered by normalization peptides)
    type: long?
  peptideEstimation__MinPeptidesPerBin:
    doc: Minimal number of peptides that are required for a bin to counted as 'covered'
    type: long?
  peptideEstimation__MinBinsFilled:
    doc: Minimal number of bins required to be covered
    type: long?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: OpenSwathRTNormalizer
doc: This tool will take a description of RT peptides and their normalized retention time to write out a transformation file on how to transform the RT space into the normalized space.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - OpenSwathRTNormalizer
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json