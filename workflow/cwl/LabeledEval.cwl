inputs:
  in:
    doc: Feature result file
    type: File
  truth:
    doc: Expected result file.
    type: File
  rt_tol:
    doc: Maximum allowed retention time deviation
    type: double?
  mz_tol:
    doc: Maximum allowed m/z deviation (divided by charge)
    type: double?
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
label: LabeledEval
doc: " Evaluation tool for isotope-labeled quantitation experiments."
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - LabeledEval
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json