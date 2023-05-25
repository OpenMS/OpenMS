inputs:
  in:
    doc: Input spectrum file
    type: File
  id:
    doc: Protein/peptide identifications file
    type: File
  out:
    doc: Output file
    type: string
  executable:
    doc: LuciPHOr2 .jar file. Provide a full or relative path, or make sure it can be found in your PATH environment.
    type: File
  fragment_method:
    doc: Fragmentation method
    type: string?
  fragment_mass_tolerance:
    doc: Tolerance of the peaks in the fragment spectrum
    type: double?
  fragment_error_units:
    doc: Unit of fragment mass tolerance
    type: string?
  min_mz:
    doc: Do not consider peaks below this value for matching fragment ions
    type: double?
  target_modifications:
    doc: List the amino acids to be searched for and their mass modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'
    type: string[]?
  neutral_losses:
    doc: "List the types of neutral losses that you want to consider. The residue field is case sensitive. For example: lower case 'sty' implies that the neutral loss can only occur if the specified modification is present. Syntax: NL = <RESDIUES> -<NEUTRAL_LOSS_MOLECULAR_FORMULA> <MASS_LOST>"
    type: string[]?
  decoy_mass:
    doc: How much to add to an amino acid to make it a decoy
    type: double?
  decoy_neutral_losses:
    doc: "For handling the neutral loss from a decoy sequence. The syntax for this is identical to that of the normal neutral losses given above except that the residue is always 'X'. Syntax: DECOY_NL = X -<NEUTRAL_LOSS_MOLECULAR_FORMULA> <MASS_LOST>"
    type: string[]?
  max_charge_state:
    doc: Do not consider PSMs with a charge state above this value
    type: long?
  max_peptide_length:
    doc: Restrict scoring to peptides with a length shorter than this value
    type: long?
  max_num_perm:
    doc: Maximum number of permutations a sequence can have
    type: long?
  modeling_score_threshold:
    doc: Minimum score a PSM needs to be considered for modeling
    type: double?
  scoring_threshold:
    doc: PSMs below this value will be discarded
    type: double?
  min_num_psms_model:
    doc: The minimum number of PSMs you need for any charge state in order to build a model for it
    type: long?
  num_threads:
    doc: For multi-threading, 0 = use all CPU found by JAVA
    type: long?
  run_mode:
    doc: "Determines how Luciphor will run: 0 = calculate FLR then rerun scoring without decoys (two iterations), 1 = Report Decoys: calculate FLR but don't rescore PSMs, all decoy hits will be reported"
    type: string?
  rt_tolerance:
    doc: Set the retention time tolerance (for the mapping of identifications to spectra in case multiple search engines were used)
    type: double?
  java_executable:
    doc: The Java executable. Usually Java is on the system PATH. If Java is not found, use this parameter to specify the full path to Java
    type: File?
  java_memory:
    doc: Maximum Java heap size (in MB)
    type: long?
  java_permgen:
    doc: Maximum Java permanent generation space (in MB); only for Java 7 and below
    type: long?
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
label: LuciphorAdapter
doc: Modification site localisation using LuciPHOr2.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - LuciphorAdapter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json