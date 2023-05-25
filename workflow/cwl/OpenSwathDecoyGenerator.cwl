inputs:
  in:
    doc: Input file
    type: File
  in_type:
    doc: "Input file type -- default: determined from file extension or content\n"
    type: string?
  out:
    doc: Output file
    type: string
  out_type:
    doc: "Output file type -- default: determined from file extension or content\n"
    type: string?
  method:
    doc: Decoy generation method
    type: string?
  decoy_tag:
    doc: decoy tag
    type: string?
  min_decoy_fraction:
    doc: Minimum fraction of decoy / target peptides and proteins
    type: double?
  aim_decoy_fraction:
    doc: Number of decoys the algorithm should generate (if unequal to 1, the algorithm will randomly select N peptides for decoy generation)
    type: double?
  shuffle_max_attempts:
    doc: "shuffle: maximum attempts to lower the amino acid sequence identity between target and decoy for the shuffle algorithm"
    type: long?
  shuffle_sequence_identity_threshold:
    doc: "shuffle: target-decoy amino acid sequence identity threshold for the shuffle algorithm"
    type: double?
  shift_precursor_mz_shift:
    doc: "shift: precursor ion MZ shift in Thomson for shift decoy method"
    type: double?
  shift_product_mz_shift:
    doc: "shift: fragment ion MZ shift in Thomson for shift decoy method"
    type: double?
  product_mz_threshold:
    doc: MZ threshold in Thomson for fragment ion annotation
    type: double?
  allowed_fragment_types:
    doc: allowed fragment types
    type: string?
  allowed_fragment_charges:
    doc: allowed fragment charge states
    type: string?
  enable_detection_specific_losses:
    doc: set this flag if specific neutral losses for detection fragment ions should be allowed
    type: boolean?
  enable_detection_unspecific_losses:
    doc: set this flag if unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) for detection fragment ions should be allowed
    type: boolean?
  switchKR:
    doc: Whether to switch terminal K and R (to achieve different precursor mass)
    type: string?
  separate:
    doc: set this flag if decoys should not be appended to targets.
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
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: OpenSwathDecoyGenerator
doc: Generates decoys according to different models for a specific TraML
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - OpenSwathDecoyGenerator
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json