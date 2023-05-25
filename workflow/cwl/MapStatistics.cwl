inputs:
  in:
    doc: Input file
    type: File
  in_type:
    doc: "Input file type -- default: determined from file extension or content"
    type: string?
  out:
    doc: Optional output txt file. If empty, the output is written to the command line.
    type: string
  n:
    doc: Report separate statistics for each of n RT slices of the map.
    type: long?
  m:
    doc: Show meta information about the whole experiment
    type: boolean?
  p:
    doc: Shows data processing information
    type: boolean?
  s:
    doc: Computes a summary statistics of intensities, qualities, and widths
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
label: MapStatistics
doc: Extract extended statistics on the features of a map for quality control.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MapStatistics
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json