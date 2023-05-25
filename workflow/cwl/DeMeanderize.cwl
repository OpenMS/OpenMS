inputs:
  in:
    doc: Input experiment file, containing the wrongly sorted spectra.
    type: File
  out:
    doc: Output experiment file with correctly sorted spectra.
    type: string
  num_spots_per_row:
    doc: Number of spots in one column, until next row is spotted.
    type: long?
  RT_distance:
    doc: RT distance between two spots which is used to calculated pseudo RT.
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
label: DeMeanderize
doc: Orders the spectra of MALDI spotting plates correctly.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - DeMeanderize
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json