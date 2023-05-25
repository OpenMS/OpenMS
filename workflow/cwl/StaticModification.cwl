inputs:
  in:
    doc: "Input: identification results"
    type: File
  out:
    doc: "Output: identification results with modifications applied"
    type: string
  mods:
    doc: List of manual modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'.
    type: string[]?
  presets:
    doc: Add predefined sets, as shortcut to manually specifying a lot of modifications.
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
label: StaticModification
doc: Applies a set of modifications to all PeptideIDs in an idXML file.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - StaticModification
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json