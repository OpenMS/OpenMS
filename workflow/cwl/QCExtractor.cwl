inputs:
  in:
    doc: Input qcml file
    type: File
  qp:
    doc: Target attachment qp.
    type: string
  run:
    doc: The file that defined the run under which the qp for the attachment is aggregated as mzML file. The file is only used to extract the run name from the file name.
    type: File?
  name:
    doc: If no file for the run was given (or if the target qp is contained in a set), at least a name of the target run/set containing the the qp for the attachment has to be given.
    type: string?
  out_csv:
    doc: Output csv formatted table.
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
  out_csv:
    type: File
    outputBinding:
      glob: $(inputs.out_csv)
label: QCExtractor
doc: Extracts a table attachment to a given qc parameter.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - QCExtractor
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json