inputs:
  in:
    doc: input files separated by blanks
    type: File[]
  out:
    doc: Output file
    type: string
  design:
    doc: input file containing the experimental design
    type: File?
  keep_subelements:
    doc: "For consensusXML input only: If set, the sub-features of the inputs are transferred to the output."
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
  algorithm__mz_unit:
    doc: Unit of m/z tolerance
    type: string?
  algorithm__nr_partitions:
    doc: Number of partitions in m/z space
    type: long?
  algorithm__warp__enabled:
    doc: Whether or not to internally warp feature RTs using LOWESS transformation before linking (reported RTs in results will always be the original RTs)
    type: string?
  algorithm__warp__rt_tol:
    doc: Width of RT tolerance window (sec)
    type: double?
  algorithm__warp__mz_tol:
    doc: m/z tolerance (in ppm or Da)
    type: double?
  algorithm__warp__max_pairwise_log_fc:
    doc: "Maximum absolute log10 fold change between two compatible signals during compatibility graph construction. Two signals from different maps will not be connected by an edge in the compatibility graph if absolute log fold change exceeds this limit (they might still end up in the same connected component, however). Note: this does not limit fold changes in the linking stage, only during RT alignment, where we try to find high-quality alignment anchor points. Setting this to a value < 0 disables the FC check."
    type: double?
  algorithm__warp__min_rel_cc_size:
    doc: Only connected components containing compatible features from at least max(2, (warp_min_occur * number_of_input_maps)) input maps are considered for computing the warping function
    type: double?
  algorithm__warp__max_nr_conflicts:
    doc: Allow up to this many conflicts (features from the same map) per connected component to be used for alignment (-1 means allow any number of conflicts)
    type: long?
  algorithm__link__rt_tol:
    doc: Width of RT tolerance window (sec)
    type: double?
  algorithm__link__mz_tol:
    doc: m/z tolerance (in ppm or Da)
    type: double?
  algorithm__link__charge_merging:
    doc: whether to disallow charge mismatches (Identical), allow to link charge zero (i.e., unknown charge state) with every charge state, or disregard charges (Any).
    type: string?
  algorithm__link__adduct_merging:
    doc: whether to only allow the same adduct for linking (Identical), also allow linking features with adduct-free ones, or disregard adducts (Any).
    type: string?
  algorithm__distance_RT__exponent:
    doc: Normalized RT differences ([0-1], relative to 'max_difference') are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)
    type: double?
  algorithm__distance_RT__weight:
    doc: Final RT distances are weighted by this factor
    type: double?
  algorithm__distance_MZ__exponent:
    doc: Normalized ([0-1], relative to 'max_difference') m/z differences are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)
    type: double?
  algorithm__distance_MZ__weight:
    doc: Final m/z distances are weighted by this factor
    type: double?
  algorithm__distance_intensity__exponent:
    doc: Differences in relative intensity ([0-1]) are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)
    type: double?
  algorithm__distance_intensity__weight:
    doc: Final intensity distances are weighted by this factor
    type: double?
  algorithm__distance_intensity__log_transform:
    doc: Log-transform intensities? If disabled, d = |int_f2 - int_f1| / int_max. If enabled, d = |log(int_f2 + 1) - log(int_f1 + 1)| / log(int_max + 1))
    type: string?
  algorithm__LOWESS__span:
    doc: Fraction of datapoints (f) to use for each local regression (determines the amount of smoothing). Choosing this parameter in the range .2 to .8 usually results in a good fit.
    type: double?
  algorithm__LOWESS__num_iterations:
    doc: Number of robustifying iterations for lowess fitting.
    type: long?
  algorithm__LOWESS__delta:
    doc: Nonnegative parameter which may be used to save computations (recommended value is 0.01 of the range of the input, e.g. for data ranging from 1000 seconds to 2000 seconds, it could be set to 10). Setting a negative value will automatically do this.
    type: double?
  algorithm__LOWESS__interpolation_type:
    doc: "Method to use for interpolation between datapoints computed by lowess. 'linear': Linear interpolation. 'cspline': Use the cubic spline for interpolation. 'akima': Use an akima spline for interpolation"
    type: string?
  algorithm__LOWESS__extrapolation_type:
    doc: "Method to use for extrapolation outside the data range. 'two-point-linear': Uses a line through the first and last point to extrapolate. 'four-point-linear': Uses a line through the first and second point to extrapolate in front and and a line through the last and second-to-last point in the end. 'global-linear': Uses a linear regression to fit a line through all data points and use it for interpolation."
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: FeatureLinkerUnlabeledKD
doc: Groups corresponding features from multiple maps.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FeatureLinkerUnlabeledKD
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json