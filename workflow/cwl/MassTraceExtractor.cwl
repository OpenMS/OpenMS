inputs:
  in:
    doc: input centroided mzML file
    type: File
  out:
    doc: output featureXML file with mass traces
    type: string
  out_type:
    doc: "output file type -- default: determined from file extension or content"
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
  algorithm__common__noise_threshold_int:
    doc: Intensity threshold below which peaks are regarded as noise.
    type: double?
  algorithm__common__chrom_peak_snr:
    doc: Minimum signal-to-noise a mass trace should have.
    type: double?
  algorithm__common__chrom_fwhm:
    doc: Expected chromatographic peak width (in seconds).
    type: double?
  algorithm__mtd__mass_error_ppm:
    doc: Allowed mass deviation (in ppm).
    type: double?
  algorithm__mtd__reestimate_mt_sd:
    doc: Enables dynamic re-estimation of m/z variance during mass trace collection stage.
    type: string?
  algorithm__mtd__quant_method:
    doc: Method of quantification for mass traces. For LC data 'area' is recommended, 'median' for direct injection data. 'max_height' simply uses the most intense peak in the trace.
    type: string?
  algorithm__mtd__trace_termination_criterion:
    doc: Termination criterion for the extension of mass traces. In 'outlier' mode, trace extension cancels if a predefined number of consecutive outliers are found (see trace_termination_outliers parameter). In 'sample_rate' mode, trace extension in both directions stops if ratio of found peaks versus visited spectra falls below the 'min_sample_rate' threshold.
    type: string?
  algorithm__mtd__trace_termination_outliers:
    doc: Mass trace extension in one direction cancels if this number of consecutive spectra with no detectable peaks is reached.
    type: long?
  algorithm__mtd__min_sample_rate:
    doc: Minimum fraction of scans along the mass trace that must contain a peak.
    type: double?
  algorithm__mtd__min_trace_length:
    doc: Minimum expected length of a mass trace (in seconds).
    type: double?
  algorithm__mtd__max_trace_length:
    doc: Maximum expected length of a mass trace (in seconds). Set to a negative value to disable maximal length check during mass trace detection.
    type: double?
  algorithm__epd__width_filtering:
    doc: Enable filtering of unlikely peak widths. The fixed setting filters out mass traces outside the [min_fwhm, max_fwhm] interval (set parameters accordingly!). The auto setting filters with the 5 and 95% quantiles of the peak width distribution.
    type: string?
  algorithm__epd__min_fwhm:
    doc: Minimum full-width-at-half-maximum of chromatographic peaks (in seconds). Ignored if parameter width_filtering is off or auto.
    type: double?
  algorithm__epd__max_fwhm:
    doc: Maximum full-width-at-half-maximum of chromatographic peaks (in seconds). Ignored if parameter width_filtering is off or auto.
    type: double?
  algorithm__epd__masstrace_snr_filtering:
    doc: Apply post-filtering by signal-to-noise ratio after smoothing.
    type: boolean?
  algorithm__epd__enabled:
    doc: Enables/disables the chromatographic peak detection of mass traces
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: MassTraceExtractor
doc: Detects mass traces in centroided LC-MS data.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MassTraceExtractor
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json