inputs:
  cv_files:
    doc: List of ontology files in OBO format.
    type: File[]
  cv_names:
    doc: List of identifiers (one for each ontology file).
    type: string[]
  mapping_file:
    doc: Mapping file in CVMapping (XML) format.
    type: File
  ignore_cv:
    doc: A list of CV identifiers which should be ignored.
    type: string[]?
  html:
    doc: Writes an HTML version of the mapping file with annotated CV terms
    type: string
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
  html:
    type: File
    outputBinding:
      glob: $(inputs.html)
label: CVInspector
doc: A tool for visualization and validation of PSI mapping and CV files.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - CVInspector
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json