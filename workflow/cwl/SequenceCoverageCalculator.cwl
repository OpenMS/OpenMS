inputs:
  in_database:
    doc: input file containing the database in FASTA format
    type: File
  in_peptides:
    doc: input file containing the identified peptides
    type: File
  out:
    doc: Optional text output file. If left out, the output is written to the command line.
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
label: SequenceCoverageCalculator
doc: Prints information about idXML files.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - SequenceCoverageCalculator
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json