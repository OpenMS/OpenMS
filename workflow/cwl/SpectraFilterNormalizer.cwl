inputs:
  in:
    doc: input file
    type: File
  out:
    doc: output file
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
  algorithm__method:
    doc: Normalize via dividing by TIC ('to_TIC') per spectrum or normalize to max. intensity of one ('to_one') per spectrum.
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: SpectraFilterNormalizer
doc: Normalizes intensity of peak spectra.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - SpectraFilterNormalizer
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json