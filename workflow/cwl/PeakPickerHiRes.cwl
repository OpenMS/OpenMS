inputs:
  in:
    doc: "input profile data file "
    type: File
  out:
    doc: "output peak file "
    type: string
  processOption:
    doc: Whether to load all data and process them in-memory or whether to process the data on the fly (lowmemory) without loading the whole file into memory first
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
  algorithm__signal_to_noise:
    doc: Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)
    type: double?
  algorithm__spacing_difference_gap:
    doc: The extension of a peak is stopped if the spacing between two subsequent data points exceeds 'spacing_difference_gap * min_spacing'. 'min_spacing' is the smaller of the two spacings from the peak apex to its two neighboring points. '0' to disable the constraint. Not applicable to chromatograms.
    type: double?
  algorithm__spacing_difference:
    doc: Maximum allowed difference between points during peak extension, in multiples of the minimal difference between the peak apex and its two neighboring points. If this difference is exceeded a missing point is assumed (see parameter 'missing'). A higher value implies a less stringent peak definition, since individual signals within the peak are allowed to be further apart. '0' to disable the constraint. Not applicable to chromatograms.
    type: double?
  algorithm__missing:
    doc: Maximum number of missing points allowed when extending a peak to the left or to the right. A missing data point occurs if the spacing between two subsequent data points exceeds 'spacing_difference * min_spacing'. 'min_spacing' is the smaller of the two spacings from the peak apex to its two neighboring points. Not applicable to chromatograms.
    type: long?
  algorithm__ms_levels:
    doc: List of MS levels for which the peak picking is applied. If empty, auto mode is enabled, all peaks which aren't picked yet will get picked. Other scans are copied to the output without changes.
    type: long[]?
  algorithm__report_FWHM:
    doc: Add metadata for FWHM (as floatDataArray named 'FWHM' or 'FWHM_ppm', depending on param 'report_FWHM_unit') for each picked peak.
    type: boolean?
  algorithm__report_FWHM_unit:
    doc: Unit of FWHM. Either absolute in the unit of input, e.g. 'm/z' for spectra, or relative as ppm (only sensible for spectra, not chromatograms).
    type: string?
  algorithm__SignalToNoise__max_intensity:
    doc: maximal intensity considered for histogram construction. By default, it will be calculated automatically (see auto_mode). Only provide this parameter if you know what you are doing (and change 'auto_mode' to '-1')! All intensities EQUAL/ABOVE 'max_intensity' will be added to the LAST histogram bin. If you choose 'max_intensity' too small, the noise estimate might be too small as well.  If chosen too big, the bins become quite large (which you could counter by increasing 'bin_count', which increases runtime). In general, the Median-S/N estimator is more robust to a manual max_intensity than the MeanIterative-S/N.
    type: long?
  algorithm__SignalToNoise__auto_max_stdev_factor:
    doc: "parameter for 'max_intensity' estimation (if 'auto_mode' == 0): mean + 'auto_max_stdev_factor' * stdev"
    type: double?
  algorithm__SignalToNoise__auto_max_percentile:
    doc: "parameter for 'max_intensity' estimation (if 'auto_mode' == 1): auto_max_percentile th percentile"
    type: long?
  algorithm__SignalToNoise__auto_mode:
    doc: "method to use to determine maximal intensity: -1 --> use 'max_intensity'; 0 --> 'auto_max_stdev_factor' method (default); 1 --> 'auto_max_percentile' method"
    type: long?
  algorithm__SignalToNoise__win_len:
    doc: window length in Thomson
    type: double?
  algorithm__SignalToNoise__bin_count:
    doc: number of bins for intensity values
    type: long?
  algorithm__SignalToNoise__min_required_elements:
    doc: minimum number of elements required in a window (otherwise it is considered sparse)
    type: long?
  algorithm__SignalToNoise__noise_for_empty_window:
    doc: noise value used for sparse windows
    type: double?
  algorithm__SignalToNoise__write_log_messages:
    doc: Write out log messages in case of sparse windows or median in rightmost histogram bin
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: PeakPickerHiRes
doc: Finds mass spectrometric peaks in profile mass spectra.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - PeakPickerHiRes
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json