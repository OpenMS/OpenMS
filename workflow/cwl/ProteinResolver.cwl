inputs:
  fasta:
    doc: Input database file
    type: File
  in:
    doc: Input file(s) holding experimental data
    type: File[]?
  in_path:
    doc: Path to idXMLs or consensusXMLs files. Ignored if 'in' is given.
    type: string?
  design:
    doc: Text file containing the experimental design. See documentation for specific format requirements
    type: File?
  protein_groups:
    doc: output file. Contains all protein groups
    type: string
  peptide_table:
    doc: output file. Contains one peptide per line and all proteins which contain that peptide
    type: string
  protein_table:
    doc: output file. Contains one protein per line
    type: string
  additional_info:
    doc: output file for additional info
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
  resolver__missed_cleavages:
    doc: Number of allowed missed cleavages
    type: long?
  resolver__min_length:
    doc: Minimum length of peptide
    type: long?
  resolver__enzyme:
    doc: Digestion enzyme
    type: string?
  designer__experiment:
    doc: Identifier for the experimental design.
    type: string?
  designer__file:
    doc: Identifier for the file name.
    type: string?
  designer__separator:
    doc: Separator, which should be used to split a row into columns
    type: string?
outputs:
  protein_groups:
    type: File
    outputBinding:
      glob: $(inputs.protein_groups)
  peptide_table:
    type: File
    outputBinding:
      glob: $(inputs.peptide_table)
  protein_table:
    type: File
    outputBinding:
      glob: $(inputs.protein_table)
  additional_info:
    type: File
    outputBinding:
      glob: $(inputs.additional_info)
label: ProteinResolver
doc: protein inference
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - ProteinResolver
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json