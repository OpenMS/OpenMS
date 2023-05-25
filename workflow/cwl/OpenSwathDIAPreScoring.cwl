inputs:
  tr:
    doc: transition file
    type: File
  swath_files:
    doc: Swath files that were used to extract the transitions. If present, SWATH specific scoring will be applied.
    type: File[]
  min_upper_edge_dist:
    doc: Minimal distance to the edge to still consider a precursor, in Thomson (only in SWATH)
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
label: OpenSwathDIAPreScoring
doc: Scoring spectra using the DIA scores.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - OpenSwathDIAPreScoring
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json