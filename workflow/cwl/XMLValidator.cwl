inputs:
  in:
    doc: file to validate
    type: File
  schema:
    doc: "schema to validate against.\nIf no schema is given, the file is validated against the latest schema of the file type."
    type: File?
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
label: XMLValidator
doc: Validates XML files against an XSD schema.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - XMLValidator
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json