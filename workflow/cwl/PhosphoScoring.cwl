inputs:
  in:
    doc: Input file with MS/MS spectra
    type: File
  id:
    doc: Identification input file which contains a search against a concatenated sequence database
    type: File
  out:
    doc: Identification output annotated with phosphorylation scores
    type: string
  fragment_mass_tolerance:
    doc: Fragment mass tolerance for spectrum comparisons
    type: double?
  fragment_mass_unit:
    doc: Unit of fragment mass tolerance
    type: string?
  max_peptide_length:
    doc: Restrict scoring to peptides with a length no greater than this value ('0' for 'no restriction')
    type: long?
  max_num_perm:
    doc: Maximum number of permutations a sequence can have to be processed ('0' for 'no restriction')
    type: long?
  unambiguous_score:
    doc: "Score to use for unambiguous assignments, where all sites on a peptide are phosphorylated. (Note: If a peptide is not phosphorylated at all, its score is set to '-1'.)"
    type: long?
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
label: PhosphoScoring
doc: Scores potential phosphorylation sites in order to localize the most probable sites.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - PhosphoScoring
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json