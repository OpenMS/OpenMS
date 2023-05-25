inputs:
  in:
    doc: Input file (data annotated with identifications)
    type: File
  out:
    doc: Output file (data without identifications). Either 'out' or 'id_out' are required. They can be used together.
    type: string
  id_out:
    doc: Output file (identifications). Either 'out' or 'id_out' are required. They can be used together.
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
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  id_out:
    type: File
    outputBinding:
      glob: $(inputs.id_out)
label: IDSplitter
doc: Splits protein/peptide identifications off of annotated data files
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - IDSplitter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json