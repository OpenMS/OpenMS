inputs:
  in:
    doc: Input file
    type: File
  sequences_file:
    doc: Fasta File
    type: File?
  out:
    doc: "Output file "
    type: string
  latex:
    doc: Indicates whether the output file format of the table should be LaTeX or TSV (default)
    type: boolean?
  p_value_dim_1:
    doc: Significance level of first dimension RT filter
    type: double?
  p_value_dim_2:
    doc: Significance level of second dimension RT filter
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
label: RTEvaluation
doc: Application that evaluates TPs (true positives), TNs, FPs, and FNs for an idXML file with predicted RTs.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - RTEvaluation
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json