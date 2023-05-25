inputs:
  in:
    doc: Input FASTA file, containing a database.
    type: File
  id:
    doc: Input file containing identified peptides and proteins.
    type: File
  method:
    doc: Switch between white-/blacklisting
    type: string?
  out:
    doc: Output FASTA file where the reduced database will be written to.
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
label: DatabaseFilter
doc: Filters a protein database (FASTA format) based on identified proteins
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - DatabaseFilter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json