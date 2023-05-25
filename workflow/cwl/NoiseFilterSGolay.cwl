inputs:
  in:
    doc: "input raw data file "
    type: File
  out:
    doc: "output raw data file "
    type: string
  processOption:
    doc: Whether to load all data and process them in-memory or whether to process the data on the fly (lowmemory) without loading the whole file into memory first
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
  algorithm__frame_length:
    doc: "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added."
    type: long?
  algorithm__polynomial_order:
    doc: Order or the polynomial that is fitted.
    type: long?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: NoiseFilterSGolay
doc: Removes noise from profile spectra by using a Savitzky Golay filter. Requires uniform (equidistant) data.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - NoiseFilterSGolay
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json