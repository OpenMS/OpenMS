inputs:
  in:
    doc: Input file (any xml file)
    type: File
  mapping_file:
    doc: Mapping file which is used to semantically validate the given XML file against this mapping file (see 'share/OpenMS/MAPPING' for templates).
    type: File
  cv:
    doc: Controlled Vocabulary files containg the CV terms (if left empty, a set of default files are used)
    type: File[]?
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
  {}
label: SemanticValidator
doc: SemanticValidator for semantically validating certain XML files.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - SemanticValidator
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json