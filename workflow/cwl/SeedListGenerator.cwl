inputs:
  in:
    doc: Input file (see below for details)
    type: File
  out_prefix:
    doc: Output file prefix
    type: string
  use_peptide_mass:
    doc: "[idXML input only] Use the monoisotopic mass of the best peptide hit for the m/z position (default: use precursor m/z)"
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
  out_prefix:
    type: File
    outputBinding:
      glob: $(inputs.out_prefix)*
label: SeedListGenerator
doc: Generates seed lists for feature detection.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - SeedListGenerator
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json