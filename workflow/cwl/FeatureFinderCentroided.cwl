inputs:
  in:
    doc: input file
    type: File
  out:
    doc: output file
    type: string
  seeds:
    doc: User specified seed list
    type: File?
  out_mzq:
    doc: Optional output file of MzQuantML.
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
  algorithm__debug:
    doc: When debug mode is activated, several files with intermediate results are written to the folder 'debug' (do not use in parallel mode).
    type: boolean?
  algorithm__intensity__bins:
    doc: "Number of bins per dimension (RT and m/z). The higher this value, the more local the intensity significance score is.\nThis parameter should be decreased, if the algorithm is used on small regions of a map."
    type: long?
  algorithm__mass_trace__mz_tolerance:
    doc: "Tolerated m/z deviation of peaks belonging to the same mass trace.\nIt should be larger than the m/z resolution of the instrument.\nThis value must be smaller than that 1/charge_high!"
    type: double?
  algorithm__mass_trace__min_spectra:
    doc: Number of spectra that have to show a similar peak mass in a mass trace.
    type: long?
  algorithm__mass_trace__max_missing:
    doc: "Number of consecutive spectra where a high mass deviation or missing peak is acceptable.\nThis parameter should be well below 'min_spectra'!"
    type: long?
  algorithm__mass_trace__slope_bound:
    doc: "The maximum slope of mass trace intensities when extending from the highest peak.\nThis parameter is important to separate overlapping elution peaks.\nIt should be increased if feature elution profiles fluctuate a lot."
    type: double?
  algorithm__isotopic_pattern__charge_low:
    doc: Lowest charge to search for.
    type: long?
  algorithm__isotopic_pattern__charge_high:
    doc: Highest charge to search for.
    type: long?
  algorithm__isotopic_pattern__mz_tolerance:
    doc: "Tolerated m/z deviation from the theoretical isotopic pattern.\nIt should be larger than the m/z resolution of the instrument.\nThis value must be smaller than that 1/charge_high!"
    type: double?
  algorithm__isotopic_pattern__intensity_percentage:
    doc: Isotopic peaks that contribute more than this percentage to the overall isotope pattern intensity must be present.
    type: double?
  algorithm__isotopic_pattern__intensity_percentage_optional:
    doc: Isotopic peaks that contribute more than this percentage to the overall isotope pattern intensity can be missing.
    type: double?
  algorithm__isotopic_pattern__optional_fit_improvement:
    doc: Minimal percental improvement of isotope fit to allow leaving out an optional peak.
    type: double?
  algorithm__isotopic_pattern__mass_window_width:
    doc: Window width in Dalton for precalculation of estimated isotope distributions.
    type: double?
  algorithm__isotopic_pattern__abundance_12C:
    doc: Rel. abundance of the light carbon. Modify if labeled.
    type: double?
  algorithm__isotopic_pattern__abundance_14N:
    doc: Rel. abundance of the light nitrogen. Modify if labeled.
    type: double?
  algorithm__seed__min_score:
    doc: "Minimum seed score a peak has to reach to be used as seed.\nThe seed score is the geometric mean of intensity score, mass trace score and isotope pattern score.\nIf your features show a large deviation from the averagene isotope distribution or from an gaussian elution profile, lower this score."
    type: double?
  algorithm__fit__max_iterations:
    doc: Maximum number of iterations of the fit.
    type: long?
  algorithm__feature__min_score:
    doc: "Feature score threshold for a feature to be reported.\nThe feature score is the geometric mean of the average relative deviation and the correlation between the model and the observed peaks."
    type: double?
  algorithm__feature__min_isotope_fit:
    doc: Minimum isotope fit of the feature before model fitting.
    type: double?
  algorithm__feature__min_trace_score:
    doc: "Trace score threshold.\nTraces below this threshold are removed after the model fitting.\nThis parameter is important for features that overlap in m/z dimension."
    type: double?
  algorithm__feature__min_rt_span:
    doc: Minimum RT span in relation to extended area that has to remain after model fitting.
    type: double?
  algorithm__feature__max_rt_span:
    doc: Maximum RT span in relation to extended area that the model is allowed to have.
    type: double?
  algorithm__feature__rt_shape:
    doc: Choose model used for RT profile fitting. If set to symmetric a gauss shape is used, in case of asymmetric an EGH shape is used.
    type: string?
  algorithm__feature__max_intersection:
    doc: Maximum allowed intersection of features.
    type: double?
  algorithm__feature__reported_mz:
    doc: "The mass type that is reported for features.\n'maximum' returns the m/z value of the highest mass trace.\n'average' returns the intensity-weighted average m/z value of all contained peaks.\n'monoisotopic' returns the monoisotopic m/z value derived from the fitted isotope model."
    type: string?
  algorithm__user-seed__rt_tolerance:
    doc: Allowed RT deviation of seeds from the user-specified seed position.
    type: double?
  algorithm__user-seed__mz_tolerance:
    doc: Allowed m/z deviation of seeds from the user-specified seed position.
    type: double?
  algorithm__user-seed__min_score:
    doc: Overwrites 'seed:min_score' for user-specified seeds. The cutoff is typically a bit lower in this case.
    type: double?
  algorithm__debug__pseudo_rt_shift:
    doc: Pseudo RT shift used when .
    type: double?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_mzq:
    type: File
    outputBinding:
      glob: $(inputs.out_mzq)
label: FeatureFinderCentroided
doc: Detects two-dimensional features in LC-MS data.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FeatureFinderCentroided
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json