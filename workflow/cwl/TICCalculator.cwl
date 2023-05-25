inputs:
  in:
    doc: Input file to convert.
    type: File
  in_type:
    doc: "Input file type -- default: determined from file extension or content\n"
    type: string?
  read_method:
    doc: Method to read the file
    type: string?
  loadData:
    doc: Whether to actually load and decode the binary data (or whether to skip decoding the binary data)
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
  {}
label: TICCalculator
doc: Calculates the TIC from a mass spectrometric raw file (useful for benchmarking).
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - TICCalculator
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json