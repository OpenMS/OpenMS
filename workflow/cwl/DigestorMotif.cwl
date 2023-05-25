inputs:
  in:
    doc: FASTA input file
    type: File
  out:
    doc: "output file (peptides)\n"
    type: string
  missed_cleavages:
    doc: the number of allowed missed cleavages
    type: long?
  mass_accuracy:
    doc: give your mass accuracy in ppb
    type: long?
  min_length:
    doc: minimum length of peptide
    type: long?
  out_option:
    doc: indicate 1 (peptide table only), 2 (statistics only) or (both peptide table + statistics)
    type: long?
  enzyme:
    doc: The enzyme used for peptide digestion.
    type: string?
  motif:
    doc: the motif for the restricted peptidome
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
label: DigestorMotif
doc: digests a protein database in-silico
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - DigestorMotif
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json