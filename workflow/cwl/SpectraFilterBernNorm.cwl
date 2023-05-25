inputs:
  in:
    doc: "input file "
    type: File
  out:
    doc: "output file "
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
  algorithm__C1:
    doc: C1 value of the normalization.
    type: double?
  algorithm__C2:
    doc: C2 value of the normalization.
    type: double?
  algorithm__threshold:
    doc: Threshold of the Bern et al. normalization.
    type: double?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: SpectraFilterBernNorm
doc: Applies thresholdfilter to peak spectra.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - SpectraFilterBernNorm
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json