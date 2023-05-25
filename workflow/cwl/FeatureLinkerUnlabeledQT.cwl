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
  algorithm__use_identifications:
    doc: Never link features that are annotated with different peptides (only the best hit per peptide identification is taken into account).
    type: boolean?
  algorithm__nr_partitions:
    doc: How many partitions in m/z space should be used for the algorithm (more partitions means faster runtime and more memory efficient execution).
    type: long?
  algorithm__min_nr_diffs_per_bin:
    doc: "If IDs are used: How many differences from matching IDs should be used to calculate a linking tolerance for unIDed features in an RT region. RT regions will be extended until that number is reached."
    type: long?
  algorithm__min_IDscore_forTolCalc:
    doc: "If IDs are used: What is the minimum score of an ID to assume a reliable match for tolerance calculation. Check your current score type!"
    type: double?
  algorithm__noID_penalty:
    doc: "If IDs are used: For the normalized distances, how high should the penalty for missing IDs be? 0 = no bias, 1 = IDs inside the max tolerances always preferred (even if much further away)."
    type: double?
  algorithm__ignore_charge:
    doc: "false [default]: pairing requires equal charge state (or at least one unknown charge '0'); true: Pairing irrespective of charge state"
    type: boolean?
  algorithm__ignore_adduct:
    doc: "true [default]: pairing requires equal adducts (or at least one without adduct annotation); true: Pairing irrespective of adducts"
    type: string?
  algorithm__distance_RT__max_difference:
    doc: Never pair features with a larger RT distance (in seconds).
    type: double?
  algorithm__distance_RT__exponent:
    doc: Normalized RT differences ([0-1], relative to 'max_difference') are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)
    type: double?
  algorithm__distance_RT__weight:
    doc: Final RT distances are weighted by this factor
    type: double?
  algorithm__distance_MZ__max_difference:
    doc: Never pair features with larger m/z distance (unit defined by 'unit')
    type: double?
  algorithm__distance_MZ__unit:
    doc: Unit of the 'max_difference' parameter
    type: string?
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
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: FeatureLinkerUnlabeledQT
doc: Groups corresponding features from multiple maps.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FeatureLinkerUnlabeledQT
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json