inputs:
  in:
    doc: "input raw data file "
    type: File
  out:
    doc: "output raw data file "
    type: string
  struc_elem_length:
    doc: Length of the structuring element (should be wider than maximal peak width - see documentation).
    type: double?
  struc_elem_unit:
    doc: Unit of 'struc_elem_length' parameter.
    type: string?
  method:
    doc: The name of the morphological filter to be applied. If you are unsure, use the default.
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
label: BaselineFilter
doc: Removes the baseline from profile spectra using a top-hat filter.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - BaselineFilter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json