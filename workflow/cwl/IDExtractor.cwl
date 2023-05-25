inputs:
  in:
    doc: input file
    type: File
  out:
    doc: output file
    type: string
  number_of_peptides:
    doc: Number of randomly chosen peptides
    type: long?
  number_of_rand_invokations:
    doc: Number of rand invocations before random draw (basically a seed)
    type: long?
  best_hits:
    doc: If this flag is set the best n peptides are chosen.
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
      glob: $(inputs.out)
label: IDExtractor
doc: Extracts 'n' peptides randomly or best 'n' from idXML files.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - IDExtractor
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json