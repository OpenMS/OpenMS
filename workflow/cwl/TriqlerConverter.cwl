inputs:
  in:
    doc: Input consensusXML with peptide intensities
    type: File
  in_design:
    doc: Experimental Design file
    type: File
  Triqler_condition:
    doc: Which column in the condition table should be used for Triqler 'Condition'
    type: string?
  reannotate_filenames:
    doc: Overwrite MS file names in consensusXML
    type: File[]?
  out:
    doc: Input CSV file for Triqler.
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
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: TriqlerConverter
doc: Converter to input for Triqler
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - TriqlerConverter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json