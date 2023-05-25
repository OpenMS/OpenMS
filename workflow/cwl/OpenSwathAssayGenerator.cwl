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
  min_transitions:
    doc: minimal number of transitions
    type: long?
  max_transitions:
    doc: maximal number of transitions
    type: long?
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
  precursor_mz_threshold:
    doc: MZ threshold in Thomson for precursor ion selection
    type: double?
  precursor_lower_mz_limit:
    doc: lower MZ limit for precursor ions
    type: double?
  precursor_upper_mz_limit:
    doc: upper MZ limit for precursor ions
    type: double?
  product_mz_threshold:
    doc: MZ threshold in Thomson for fragment ion annotation
    type: double?
  product_lower_mz_limit:
    doc: lower MZ limit for fragment ions
    type: double?
  product_upper_mz_limit:
    doc: upper MZ limit for fragment ions
    type: double?
  swath_windows_file:
    doc: "Tab separated file containing the SWATH windows for exclusion of fragment ions falling into the precursor isolation window: lower_offset upper_offset \\newline 400 425 \\newline ... Note that the first line is a header and will be skipped."
    type: File?
  unimod_file:
    doc: (Modified) Unimod XML file (http://www.unimod.org/xml/unimod.xml) describing residue modifiability
    type: File?
  enable_ipf:
    doc: "IPF: set this flag if identification transitions should be generated for IPF. Note: Requires setting 'unimod_file'."
    type: boolean?
  max_num_alternative_localizations:
    doc: "IPF: maximum number of site-localization permutations"
    type: long?
  disable_identification_ms2_precursors:
    doc: "IPF: set this flag if MS2-level precursor ions for identification should not be allowed for extraction of the precursor signal from the fragment ion data (MS2-level)."
    type: boolean?
  disable_identification_specific_losses:
    doc: "IPF: set this flag if specific neutral losses for identification fragment ions should not be allowed"
    type: boolean?
  enable_identification_unspecific_losses:
    doc: "IPF: set this flag if unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) for identification fragment ions should be allowed"
    type: boolean?
  enable_swath_specifity:
    doc: "IPF: set this flag if identification transitions without precursor specificity (i.e. across whole precursor isolation window instead of precursor MZ) should be generated."
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
label: OpenSwathAssayGenerator
doc: Generates assays according to different models for a specific TraML
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - OpenSwathAssayGenerator
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json