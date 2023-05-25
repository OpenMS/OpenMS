inputs:
  in:
    doc: Input file
    type: File
  out:
    doc: Output file
    type: string
  new_score:
    doc: Name of the meta value to use as the new score
    type: string?
  new_score_orientation:
    doc: Orientation of the new score (are higher or lower values better?)
    type: string?
  new_score_type:
    doc: "Name to use as the type of the new score (default: same as 'new_score')"
    type: string?
  old_score:
    doc: "Name to use for the meta value storing the old score (default: old score type)"
    type: string?
  proteins:
    doc: Apply to protein scores instead of PSM scores
    type: boolean?
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
label: IDScoreSwitcher
doc: Switches between different scores of peptide or protein hits in identification data
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - IDScoreSwitcher
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json