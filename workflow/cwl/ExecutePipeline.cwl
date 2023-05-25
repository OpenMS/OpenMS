inputs:
  in:
    doc: The workflow to be executed.
    type: File
  out_dir:
    doc: "Directory for output files (default: user's home directory)"
    type: string?
  resource_file:
    doc: A TOPPAS resource file (*.trf) specifying the files this workflow is to be applied to
    type: string?
  num_jobs:
    doc: Maximum number of jobs running in parallel
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
outputs:
  {}
label: ExecutePipeline
doc: Executes workflows created by TOPPAS.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - ExecutePipeline
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json