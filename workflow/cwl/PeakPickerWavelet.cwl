inputs:
  in:
    doc: "input profile data file "
    type: File
  out:
    doc: "output peak file "
    type: string
  write_peak_meta_data:
    doc: "Write additional information about the picked peaks (maximal intensity, left and right area...) into the mzML-file. Attention: this can blow up files, since seven arrays are stored per spectrum!"
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
  algorithm__signal_to_noise:
    doc: Minimal signal to noise ratio for a peak to be picked.
    type: double?
  algorithm__centroid_percentage:
    doc: Percentage of the maximum height that the raw data points must exceed to be taken into account for the calculation of the centroid. If it is 1 the centroid position corresponds to the position of the highest intensity.
    type: double?
  algorithm__peak_width:
    doc: Approximate fwhm of the peaks.
    type: double?
  algorithm__estimate_peak_width:
    doc: "Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored."
    type: boolean?
  algorithm__fwhm_lower_bound_factor:
    doc: Factor that calculates the minimal fwhm value from the peak_width. All peaks with width smaller than fwhm_bound_factor * peak_width are discarded.
    type: double?
  algorithm__fwhm_upper_bound_factor:
    doc: Factor that calculates the maximal fwhm value from the peak_width. All peaks with width greater than fwhm_upper_bound_factor * peak_width are discarded.
    type: double?
  algorithm__optimization:
    doc: If the peak parameters position, intensity and left/right widthshall be optimized set optimization to one_dimensional or two_dimensional.
    type: string?
  algorithm__thresholds__peak_bound:
    doc: Minimal peak intensity.
    type: double?
  algorithm__thresholds__peak_bound_ms2_level:
    doc: Minimal peak intensity for MS/MS peaks.
    type: double?
  algorithm__thresholds__correlation:
    doc: minimal correlation of a peak and the raw signal. If a peak has a lower correlation it is skipped.
    type: double?
  algorithm__thresholds__noise_level:
    doc: noise level for the search of the peak endpoints.
    type: double?
  algorithm__thresholds__search_radius:
    doc: search radius for the search of the maximum in the signal after a maximum in the cwt was found
    type: long?
  algorithm__wavelet_transform__spacing:
    doc: Spacing of the CWT. Note that the accuracy of the picked peak's centroid position depends in the Raw data spacing, i.e., 50% of raw peak distance at most.
    type: double?
  algorithm__optimization__iterations:
    doc: maximal number of iterations for the fitting step
    type: long?
  algorithm__optimization__penalties__position:
    doc: "penalty term for the fitting of the position:If it differs too much from the initial one it can be penalized "
    type: double?
  algorithm__optimization__penalties__left_width:
    doc: penalty term for the fitting of the left width:If the left width differs too much from the initial one during the fitting it can be penalized.
    type: double?
  algorithm__optimization__penalties__right_width:
    doc: penalty term for the fitting of the right width:If the right width differs too much from the initial one during the fitting it can be penalized.
    type: double?
  algorithm__optimization__penalties__height:
    doc: penalty term for the fitting of the intensity (only used in 2D Optimization):If it gets negative during the fitting it can be penalized.
    type: double?
  algorithm__optimization__2d__tolerance_mz:
    doc: mz tolerance for cluster construction
    type: double?
  algorithm__optimization__2d__max_peak_distance:
    doc: maximal peak distance in mz in a cluster
    type: double?
  algorithm__deconvolution__deconvolution:
    doc: If you want heavily overlapping peaks to be separated set this value to "true"
    type: boolean?
  algorithm__deconvolution__asym_threshold:
    doc: If the symmetry of a peak is smaller than asym_thresholds it is assumed that it consists of more than one peak and the deconvolution procedure is started.
    type: double?
  algorithm__deconvolution__left_width:
    doc: 1/left_width is the initial value for the left width of the peaks found in the deconvolution step.
    type: double?
  algorithm__deconvolution__right_width:
    doc: 1/right_width is the initial value for the right width of the peaks found in the deconvolution step.
    type: double?
  algorithm__deconvolution__scaling:
    doc: Initial scaling of the cwt used in the separation of heavily overlapping peaks. The initial value is used for charge 1, for higher charges it is adapted to scaling/charge.
    type: double?
  algorithm__deconvolution__fitting__fwhm_threshold:
    doc: If the FWHM of a peak is higher than 'fwhm_thresholds' it is assumed that it consists of more than one peak and the deconvolution procedure is started.
    type: double?
  algorithm__deconvolution__fitting__eps_abs:
    doc: if the absolute error gets smaller than this value the fitting is stopped.
    type: double?
  algorithm__deconvolution__fitting__eps_rel:
    doc: if the relative error gets smaller than this value the fitting is stopped.
    type: double?
  algorithm__deconvolution__fitting__max_iteration:
    doc: maximal number of iterations for the fitting step
    type: long?
  algorithm__deconvolution__fitting__penalties__position:
    doc: penalty term for the fitting of the peak position:If the position changes more than 0.5Da during the fitting it can be penalized as well as discrepancies of the peptide mass rule.
    type: double?
  algorithm__deconvolution__fitting__penalties__height:
    doc: penalty term for the fitting of the intensity:If it gets negative during the fitting it can be penalized.
    type: double?
  algorithm__deconvolution__fitting__penalties__left_width:
    doc: penalty term for the fitting of the left width:If the left width gets too broad or negative during the fitting it can be penalized.
    type: double?
  algorithm__deconvolution__fitting__penalties__right_width:
    doc: penalty term for the fitting of the right width:If the right width gets too broad or negative during the fitting it can be penalized.
    type: double?
  algorithm__SignalToNoiseEstimationParameter__max_intensity:
    doc: maximal intensity considered for histogram construction. By default, it will be calculated automatically (see auto_mode). Only provide this parameter if you know what you are doing (and change 'auto_mode' to '-1')! All intensities EQUAL/ABOVE 'max_intensity' will not be added to the histogram. If you choose 'max_intensity' too small, the noise estimate might be too small as well. If chosen too big, the bins become quite large (which you could counter by increasing 'bin_count', which increases runtime).
    type: long?
  algorithm__SignalToNoiseEstimationParameter__auto_max_stdev_factor:
    doc: "parameter for 'max_intensity' estimation (if 'auto_mode' == 0): mean + 'auto_max_stdev_factor' * stdev"
    type: double?
  algorithm__SignalToNoiseEstimationParameter__auto_max_percentile:
    doc: "parameter for 'max_intensity' estimation (if 'auto_mode' == 1): auto_max_percentile th percentile"
    type: long?
  algorithm__SignalToNoiseEstimationParameter__auto_mode:
    doc: "method to use to determine maximal intensity: -1 --> use 'max_intensity'; 0 --> 'auto_max_stdev_factor' method (default); 1 --> 'auto_max_percentile' method"
    type: long?
  algorithm__SignalToNoiseEstimationParameter__win_len:
    doc: window length in Thomson
    type: double?
  algorithm__SignalToNoiseEstimationParameter__bin_count:
    doc: number of bins for intensity values
    type: long?
  algorithm__SignalToNoiseEstimationParameter__stdev_mp:
    doc: multiplier for stdev
    type: double?
  algorithm__SignalToNoiseEstimationParameter__min_required_elements:
    doc: minimum number of elements required in a window (otherwise it is considered sparse)
    type: long?
  algorithm__SignalToNoiseEstimationParameter__noise_for_empty_window:
    doc: noise value used for sparse windows
    type: double?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: PeakPickerWavelet
doc: Finds mass spectrometric peaks in profile mass spectra.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - PeakPickerWavelet
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json