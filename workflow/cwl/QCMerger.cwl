inputs:
  in:
    doc: List of qcml files to be merged.
    type: File[]
  out:
    doc: Output extended/reduced qcML file
    type: string
  setname:
    doc: Use only when all given qcml files belong to one set, which will be held under the given name.
    type: string?
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
label: QCMerger
doc: Merges two qcml files together.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - QCMerger
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json