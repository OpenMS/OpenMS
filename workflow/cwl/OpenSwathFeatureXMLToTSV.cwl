inputs:
  in:
    doc: Input files separated by blank
    type: File[]
  tr:
    doc: TraML transition file
    type: File
  out:
    doc: tsv output file (mProphet compatible)
    type: string
  short_format:
    doc: Whether to write short (one peptide per line) or long format (one transition per line).
    type: boolean?
  best_scoring_peptide:
    doc: If only the best scoring feature per peptide should be printed, give the variable name
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
label: OpenSwathFeatureXMLToTSV
doc: Converts a featureXML to a mProphet tsv.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - OpenSwathFeatureXMLToTSV
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json