inputs:
  in:
    doc: Input files to align (all must have the same file type)
    type: File[]
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
  reference__file:
    doc: File to use as reference (same file format as input files required)
    type: File?
  reference__index:
    doc: "Use one of the input files as reference ('1' for the first file, etc.).\nIf '0', no explicit reference is set - the algorithm will select a reference."
    type: long?
  algorithm__max_num_peaks_considered:
    doc: The maximal number of peaks/features to be considered per map. To use all, set to '-1'.
    type: long?
  algorithm__superimposer__mz_pair_max_distance:
    doc: Maximum of m/z deviation of corresponding elements in different maps.  This condition applies to the pairs considered in hashing.
    type: double?
  algorithm__superimposer__rt_pair_distance_fraction:
    doc: "Within each of the two maps, the pairs considered for pose clustering must be separated by at least this fraction of the total elution time interval (i.e., max - min).  "
    type: double?
  algorithm__superimposer__num_used_points:
    doc: Maximum number of elements considered in each map (selected by intensity).  Use this to reduce the running time and to disregard weak signals during alignment.  For using all points, set this to -1.
    type: long?
  algorithm__superimposer__scaling_bucket_size:
    doc: The scaling of the retention time interval is being hashed into buckets of this size during pose clustering.  A good choice for this would be a bit smaller than the error you would expect from repeated runs.
    type: double?
  algorithm__superimposer__shift_bucket_size:
    doc: The shift at the lower (respectively, higher) end of the retention time interval is being hashed into buckets of this size during pose clustering.  A good choice for this would be about the time between consecutive MS scans.
    type: double?
  algorithm__superimposer__max_shift:
    doc: Maximal shift which is considered during histogramming (in seconds).  This applies for both directions.
    type: double?
  algorithm__superimposer__max_scaling:
    doc: Maximal scaling which is considered during histogramming.  The minimal scaling is the reciprocal of this.
    type: double?
  algorithm__superimposer__dump_buckets:
    doc: "[DEBUG] If non-empty, base filename where hash table buckets will be dumped to.  A serial number for each invocation will be appended automatically."
    type: string?
  algorithm__superimposer__dump_pairs:
    doc: "[DEBUG] If non-empty, base filename where the individual hashed pairs will be dumped to (large!).  A serial number for each invocation will be appended automatically."
    type: string?
  algorithm__pairfinder__second_nearest_gap:
    doc: Only link features whose distance to the second nearest neighbors (for both sides) is larger by 'second_nearest_gap' than the distance between the matched pair itself.
    type: double?
  algorithm__pairfinder__use_identifications:
    doc: Never link features that are annotated with different peptides (features without ID's always match; only the best hit per peptide identification is considered).
    type: boolean?
  algorithm__pairfinder__ignore_charge:
    doc: "false [default]: pairing requires equal charge state (or at least one unknown charge '0'); true: Pairing irrespective of charge state"
    type: boolean?
  algorithm__pairfinder__ignore_adduct:
    doc: "true [default]: pairing requires equal adducts (or at least one without adduct annotation); true: Pairing irrespective of adducts"
    type: string?
  algorithm__pairfinder__distance_RT__max_difference:
    doc: Never pair features with a larger RT distance (in seconds).
    type: double?
  algorithm__pairfinder__distance_RT__exponent:
    doc: Normalized RT differences ([0-1], relative to 'max_difference') are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)
    type: double?
  algorithm__pairfinder__distance_RT__weight:
    doc: Final RT distances are weighted by this factor
    type: double?
  algorithm__pairfinder__distance_MZ__max_difference:
    doc: Never pair features with larger m/z distance (unit defined by 'unit')
    type: double?
  algorithm__pairfinder__distance_MZ__unit:
    doc: Unit of the 'max_difference' parameter
    type: string?
  algorithm__pairfinder__distance_MZ__exponent:
    doc: Normalized ([0-1], relative to 'max_difference') m/z differences are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)
    type: double?
  algorithm__pairfinder__distance_MZ__weight:
    doc: Final m/z distances are weighted by this factor
    type: double?
  algorithm__pairfinder__distance_intensity__exponent:
    doc: Differences in relative intensity ([0-1]) are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)
    type: double?
  algorithm__pairfinder__distance_intensity__weight:
    doc: Final intensity distances are weighted by this factor
    type: double?
  algorithm__pairfinder__distance_intensity__log_transform:
    doc: Log-transform intensities? If disabled, d = |int_f2 - int_f1| / int_max. If enabled, d = |log(int_f2 + 1) - log(int_f1 + 1)| / log(int_max + 1))
    type: string?
outputs:
  {}
label: MapAlignerPoseClustering
doc: Corrects retention time distortions between maps using a pose clustering approach.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MapAlignerPoseClustering
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json