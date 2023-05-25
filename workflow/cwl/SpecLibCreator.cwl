inputs:
  info:
    doc: Holds id, peptide, retention time etc.
    type: File
  itemseperator:
    doc: " Separator between items. e.g. ,"
    type: string?
  itemenclosed:
    doc: "'true' or 'false' if true every item is enclosed e.g. '$peptide$,$run$..."
    type: boolean?
  spec:
    doc: spectra
    type: File
  out:
    doc: output MSP formatted spectra library
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
label: SpecLibCreator
doc: Creates an MSP formatted spectral library.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - SpecLibCreator
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json