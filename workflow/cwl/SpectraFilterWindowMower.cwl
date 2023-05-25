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
  algorithm__windowsize:
    doc: The size of the sliding window along the m/z axis.
    type: double?
  algorithm__peakcount:
    doc: The number of peaks that should be kept.
    type: long?
  algorithm__movetype:
    doc: Whether sliding window (one peak steps) or jumping window (window size steps) should be used.
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: SpectraFilterWindowMower
doc: Applies thresholdfilter to peak spectra.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - SpectraFilterWindowMower
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json