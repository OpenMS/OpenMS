inputs:
  in:
    doc: Input file containing RNA sequences
    type: File
  out:
    doc: Output file containing sequence fragments
    type: string
  missed_cleavages:
    doc: The number of allowed missed cleavages
    type: long?
  min_length:
    doc: Minimum length of a fragment
    type: long?
  max_length:
    doc: Maximum length of a fragment
    type: long?
  enzyme:
    doc: Digestion enzyme (RNase)
    type: string?
  unique:
    doc: Report each unique sequence fragment only once
    type: boolean?
  cdna:
    doc: Input file contains cDNA sequences - replace 'T' with 'U')
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
label: RNADigestor
doc: Digests an RNA sequence database in-silico.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - RNADigestor
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json