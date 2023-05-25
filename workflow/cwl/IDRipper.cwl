inputs:
  in:
    doc: Input file, in which the protein/peptide identifications must be tagged with 'file_origin'
    type: File
  out:
    doc: Path to the output directory to write the ripped files to.
    type: string
  numeric_filenames:
    doc: Do not infer output filenames from spectra_data or file_origin but use the input filename with numeric suffixes.
    type: boolean?
  split_ident_runs:
    doc: Split different identification runs into separate files.
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
      glob: $(inputs.out)*
label: IDRipper
doc: Split protein/peptide identification file into several files according to identification run and annotated file origin.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - IDRipper
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json