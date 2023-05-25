inputs:
  in:
    doc: Input file containing the transformation description
    type: File
  out:
    doc: Output file containing original and transformed values; if empty, output is written to the screen
    type: string
  min:
    doc: Minimum value to transform
    type: double?
  max:
    doc: Maximum value to transform (if at or below 'min', select a suitable maximum based on the transformation description)
    type: double?
  step:
    doc: Step size between 'min' and 'max'
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
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: TransformationEvaluation
doc: Applies a transformation to a range of values
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - TransformationEvaluation
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json