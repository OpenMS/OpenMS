inputs:
  in:
    doc: INI/TOPPAS files that need updating.
    type: File[]
  i:
    doc: "in-place: Override given INI/TOPPAS files with new content (not compatible with -out)"
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
  {}
label: INIUpdater
doc: Update INI and TOPPAS files to new OpenMS version.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - INIUpdater
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json