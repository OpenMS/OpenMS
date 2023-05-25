inputs:
  executable:
    doc: novor.jar
    type: File?
  in:
    doc: MzML Input file
    type: File
  out:
    doc: Novor idXML output
    type: string
  enzyme:
    doc: "Digestion enzyme - currently only Trypsin is supported "
    type: string?
  fragmentation:
    doc: Fragmentation method
    type: string?
  massAnalyzer:
    doc: MassAnalyzer e.g. (Oritrap CID-Trap, CID-FT, HCD-FT; QTof CID-TOF)
    type: string?
  fragment_mass_tolerance:
    doc: Fragmentation error tolerance  (Da)
    type: double?
  precursor_mass_tolerance:
    doc: Precursor error tolerance  (ppm or Da)
    type: double?
  precursor_error_units:
    doc: Unit of precursor mass tolerance
    type: string?
  variable_modifications:
    doc: Variable modifications
    type: string[]?
  fixed_modifications:
    doc: Fixed modifications
    type: string[]?
  forbiddenResidues:
    doc: Forbidden Resiudes
    type: string[]?
  novorFile:
    doc: File to introduce customized algorithm parameters for advanced users (otional .novor file)
    type: File?
  java_executable:
    doc: The Java executable. Usually Java is on the system PATH. If Java is not found, use this parameter to specify the full path to Java
    type: File?
  java_memory:
    doc: Maximum Java heap size (in MB)
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
label: NovorAdapter
doc: Performs de novo sequencing of peptides from MS/MS data with Novor.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - NovorAdapter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json