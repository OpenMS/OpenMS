inputs:
  in:
    doc: Input file with peptide sequences and optionally charge numbers (mutually exclusive to 'in_seq')
    type: File?
  in_seq:
    doc: List of peptide sequences (mutually exclusive to 'in')
    type: string[]?
  out:
    doc: Output file; if empty, output is written to the screen
    type: string
  charge:
    doc: List of charge states; required if 'in_seq' is given
    type: long[]?
  format:
    doc: "Output format ('list': human-readable list, 'table': CSV-like table, 'mass_only': mass values only, 'mz_only': m/z values only)\n"
    type: string?
  average_mass:
    doc: Compute average (instead of monoisotopic) peptide masses
    type: boolean?
  fragment_type:
    doc: "For what type of sequence/fragment the mass should be computed\n"
    type: string?
  separator:
    doc: Field separator for 'table' output format; by default, the 'tab' character is used
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
label: MassCalculator
doc: Calculates masses and mass-to-charge ratios of peptide sequences
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MassCalculator
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json