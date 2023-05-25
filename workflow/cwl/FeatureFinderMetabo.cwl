inputs:
  in:
    doc: Centroided mzML file
    type: File
  out:
    doc: FeatureXML file with metabolite features
    type: string
  out_chrom:
    doc: Optional mzML file with chromatograms
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
  algorithm__epd__enabled:
    doc: Enable splitting of isobaric mass traces by chromatographic peak detection. Disable for direct injection.
    type: string?
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
  algorithm__ffm__local_rt_range:
    doc: RT range where to look for coeluting mass traces
    type: double?
  algorithm__ffm__local_mz_range:
    doc: MZ range where to look for isotopic mass traces
    type: double?
  algorithm__ffm__charge_lower_bound:
    doc: Lowest charge state to consider
    type: long?
  algorithm__ffm__charge_upper_bound:
    doc: Highest charge state to consider
    type: long?
  algorithm__ffm__report_summed_ints:
    doc: Set to true for a feature intensity summed up over all traces rather than using monoisotopic trace intensity alone.
    type: string?
  algorithm__ffm__enable_RT_filtering:
    doc: Require sufficient overlap in RT while assembling mass traces. Disable for direct injection data..
    type: string?
  algorithm__ffm__isotope_filtering_model:
    doc: Remove/score candidate assemblies based on isotope intensities. SVM isotope models for metabolites were trained with either 2% or 5% RMS error. For peptides, an averagine cosine scoring is used. Select the appropriate noise model according to the quality of measurement or MS device.
    type: string?
  algorithm__ffm__mz_scoring_13C:
    doc: Use the 13C isotope peak position (~1.003355 Da) as the expected shift in m/z for isotope mass traces (highly recommended for lipidomics!). Disable for general metabolites (as described in Kenar et al. 2014, MCP.).
    type: string?
  algorithm__ffm__use_smoothed_intensities:
    doc: Use LOWESS intensities instead of raw intensities.
    type: string?
  algorithm__ffm__report_convex_hulls:
    doc: Augment each reported feature with the convex hull of the underlying mass traces (increases featureXML file size considerably).
    type: string?
  algorithm__ffm__remove_single_traces:
    doc: Remove unassembled traces (single traces).
    type: string?
  algorithm__ffm__mz_scoring_by_elements:
    doc: Use the m/z range of the assumed elements to detect isotope peaks. A expected m/z range is computed from the isotopes of the assumed elements. If enabled, this ignores 'mz_scoring_13C'
    type: string?
  algorithm__ffm__elements:
    doc: Elements assumes to be present in the sample (this influences isotope detection).
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_chrom:
    type: File
    outputBinding:
      glob: $(inputs.out_chrom)
label: FeatureFinderMetabo
doc: Assembles metabolite features from centroided (LC-)MS data using the mass trace approach.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FeatureFinderMetabo
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json