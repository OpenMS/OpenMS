inputs:
  in:
    doc: Feature input file, which contains the data to be tested against the truth file.
    type: File
  truth:
    doc: Truth feature file that defines what features should be found.
    type: File
  rt_tol:
    doc: Allowed tolerance of RT relative to average feature RT span.
    type: double?
  rt_tol_abs:
    doc: Allowed absolute tolerance of RT (overwrites 'rt_tol' if set above zero).
    type: double?
  mz_tol:
    doc: Allowed tolerance in m/z (is divided by charge).
    type: double?
  out:
    doc: Feature output file. If given, an annotated input file is written.
    type: string
  abort_reasons:
    doc: Feature file containing seeds with abort reasons.
    type: File?
  out_roc:
    doc: If given, a ROC curve file is created (ROC points based on intensity threshold)
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
  out_roc:
    type: File
    outputBinding:
      glob: $(inputs.out_roc)
label: FFEval
doc: Evaluation tool for feature detection algorithms.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FFEval
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json