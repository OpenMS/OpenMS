inputs:
  in:
    doc: "input file "
    type: File
  out:
    doc: base name of DTA output files (RT, m/z and extension are appended)
    type: string
  mz:
    doc: "m/z range of precursor peaks to extract.\nThis option is ignored for MS level 1"
    type: string?
  rt:
    doc: retention time range of spectra to extract
    type: string?
  level:
    doc: MS levels to extract
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
  {}
label: DTAExtractor
doc: Extracts spectra of an MS run file to several files in DTA format.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - DTAExtractor
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json