inputs:
  in:
    doc: Input file
    type: File
  tr:
    doc: transition file ('TraML' or 'csv')
    type: File
  out:
    doc: output file
    type: string
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
  algorithm__stop_after_feature:
    doc: Stop finding after feature (ordered by intensity; -1 means do not stop).
    type: long?
  algorithm__stop_after_intensity_ratio:
    doc: Stop after reaching intensity ratio
    type: double?
  algorithm__min_peak_width:
    doc: Minimal peak width (s), discard all peaks below this value (-1 means no action).
    type: double?
  algorithm__peak_integration:
    doc: Calculate the peak area and height either the smoothed or the raw chromatogram data.
    type: string?
  algorithm__background_subtraction:
    doc: "Remove background from peak signal using estimated noise levels. The 'original' method is only provided for historical purposes, please use the 'exact' method and set parameters using the PeakIntegrator: settings. The same original or smoothed chromatogram specified by peak_integration will be used for background estimation."
    type: string?
  algorithm__recalculate_peaks:
    doc: Tries to get better peak picking by looking at peak consistency of all picked peaks. Tries to use the consensus (median) peak border if the variation within the picked peaks is too large.
    type: boolean?
  algorithm__use_precursors:
    doc: Use precursor chromatogram for peak picking (note that this may lead to precursor signal driving the peak picking)
    type: boolean?
  algorithm__use_consensus:
    doc: Use consensus peak boundaries when computing transition group picking (if false, compute independent peak boundaries for each transition)
    type: string?
  algorithm__recalculate_peaks_max_z:
    doc: Determines the maximal Z-Score (difference measured in standard deviations) that is considered too large for peak boundaries. If the Z-Score is above this value, the median is used for peak boundaries (default value 1.0).
    type: double?
  algorithm__minimal_quality:
    doc: Only if compute_peak_quality is set, this parameter will not consider peaks below this quality threshold
    type: double?
  algorithm__resample_boundary:
    doc: For computing peak quality, how many extra seconds should be sample left and right of the actual peak
    type: double?
  algorithm__compute_peak_quality:
    doc: Tries to compute a quality value for each peakgroup and detect outlier transitions. The resulting score is centered around zero and values above 0 are generally good and below -1 or -2 are usually bad.
    type: boolean?
  algorithm__compute_peak_shape_metrics:
    doc: Calculates various peak shape metrics (e.g., tailing) that can be used for downstream QC/QA.
    type: boolean?
  algorithm__compute_total_mi:
    doc: Compute mutual information metrics for individual transitions that can be used for OpenSWATH/IPF scoring.
    type: boolean?
  algorithm__boundary_selection_method:
    doc: Method to use when selecting the best boundaries for peaks.
    type: string?
  algorithm__PeakPickerMRM__sgolay_frame_length:
    doc: "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added."
    type: long?
  algorithm__PeakPickerMRM__sgolay_polynomial_order:
    doc: Order of the polynomial that is fitted.
    type: long?
  algorithm__PeakPickerMRM__gauss_width:
    doc: Gaussian width in seconds, estimated peak size.
    type: double?
  algorithm__PeakPickerMRM__use_gauss:
    doc: Use Gaussian filter for smoothing (alternative is Savitzky-Golay filter)
    type: string?
  algorithm__PeakPickerMRM__peak_width:
    doc: Force a certain minimal peak_width on the data (e.g. extend the peak at least by this amount on both sides) in seconds. -1 turns this feature off.
    type: double?
  algorithm__PeakPickerMRM__signal_to_noise:
    doc: Signal-to-noise threshold at which a peak will not be extended any more. Note that setting this too high (e.g. 1.0) can lead to peaks whose flanks are not fully captured.
    type: double?
  algorithm__PeakPickerMRM__sn_win_len:
    doc: Signal to noise window length.
    type: double?
  algorithm__PeakPickerMRM__sn_bin_count:
    doc: Signal to noise bin count.
    type: long?
  algorithm__PeakPickerMRM__write_sn_log_messages:
    doc: Write out log messages of the signal-to-noise estimator in case of sparse windows or median in rightmost histogram bin
    type: boolean?
  algorithm__PeakPickerMRM__remove_overlapping_peaks:
    doc: Try to remove overlapping peaks during peak picking
    type: string?
  algorithm__PeakPickerMRM__method:
    doc: Which method to choose for chromatographic peak-picking (OpenSWATH legacy on raw data, corrected picking on smoothed chromatogram or Crawdad on smoothed chromatogram).
    type: string?
  algorithm__PeakIntegrator__integration_type:
    doc: The integration technique to use in integratePeak() and estimateBackground() which uses either the summed intensity, integration by Simpson's rule or trapezoidal integration.
    type: string?
  algorithm__PeakIntegrator__baseline_type:
    doc: The baseline type to use in estimateBackground() based on the peak boundaries. A rectangular baseline shape is computed based either on the minimal intensity of the peak boundaries, the maximum intensity or the average intensity (base_to_base).
    type: string?
  algorithm__PeakIntegrator__fit_EMG:
    doc: Fit the chromatogram/spectrum to the EMG peak model.
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: MRMTransitionGroupPicker
doc: Picks peaks in SRM/MRM chromatograms.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MRMTransitionGroupPicker
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json