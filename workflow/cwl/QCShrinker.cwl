inputs:
  in:
    doc: Input qcml file
    type: File
  qp_accessions:
    doc: A list of cv accessions that should be removed. If empty, the usual suspects will be removed!
    type: string[]?
  name:
    doc: The name of the target run or set that contains the requested quality parameter.
    type: string?
  run:
    doc: The file from which the name of the target run that contains the requested quality parameter is taken. This overrides the name parameter!
    type: File?
  out:
    doc: Output extended/reduced qcML file
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
label: QCShrinker
doc: This application is used to remove the verbose table attachments from a qcml file that are not needed anymore, e.g. for a final report.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - QCShrinker
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json