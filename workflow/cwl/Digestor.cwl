inputs:
  in:
    doc: input file
    type: File
  out:
    doc: Output file (peptides)
    type: string
  out_type:
    doc: Set this if you cannot control the filename of 'out', e.g., in TOPPAS.
    type: string?
  missed_cleavages:
    doc: The number of allowed missed cleavages
    type: long?
  min_length:
    doc: Minimum length of peptide
    type: long?
  max_length:
    doc: Maximum length of peptide
    type: long?
  enzyme:
    doc: The type of digestion enzyme
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
  FASTA__ID:
    doc: "Identifier to use for each peptide: copy from parent protein (parent); a consecutive number (number); parent ID + consecutive number (both)"
    type: string?
  FASTA__description:
    doc: Keep or remove the (possibly lengthy) FASTA header description. Keeping it can increase resulting FASTA file significantly.
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: Digestor
doc: Digests a protein database in-silico.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - Digestor
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json