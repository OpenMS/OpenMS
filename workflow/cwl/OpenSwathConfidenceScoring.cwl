inputs:
  in:
    doc: Input file (OpenSwath results)
    type: File
  lib:
    doc: Assay library
    type: File
  out:
    doc: Output file (results with confidence scores)
    type: string
  trafo:
    doc: Retention time transformation
    type: File?
  decoys:
    doc: Number of decoy assays to select from the library for every true assay (0 for "all")
    type: long?
  transitions:
    doc: Number of transitions per feature to consider (highest intensities first; 0 for "all")
    type: long?
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
  GLM__intercept:
    doc: Intercept term
    type: double?
  GLM__delta_rt:
    doc: Coefficient of retention time difference
    type: double?
  GLM__dist_int:
    doc: Coefficient of intensity distance
    type: double?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: OpenSwathConfidenceScoring
doc: Compute confidence scores for OpenSwath results
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - OpenSwathConfidenceScoring
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json