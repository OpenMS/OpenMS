inputs:
  in:
    doc: Input idXML file containing the identifications.
    type: File
  fasta:
    doc: Input sequence database in FASTA format. Leave empty for using the same DB as used for the input idXML (this might fail). Non-existing relative filenames are looked up via 'OpenMS.ini:id_db_dir'
    type: File?
  out:
    doc: Output idXML file.
    type: string
  decoy_string:
    doc: String that was appended (or prefixed - see 'decoy_string_position' flag below) to the accessions in the protein database to indicate decoy proteins. If empty (default), it's determined automatically (checking for common terms, both as prefix and suffix).
    type: string?
  decoy_string_position:
    doc: Is the 'decoy_string' prepended (prefix) or appended (suffix) to the protein accession? (ignored if decoy_string is empty)
    type: string?
  missing_decoy_action:
    doc: "Action to take if NO peptide was assigned to a decoy protein (which indicates wrong database or decoy string): 'error' (exit with error, no output), 'warn' (exit with success, warning message), 'silent' (no action is taken, not even a warning)"
    type: string?
  write_protein_sequence:
    doc: If set, the protein sequences are stored as well.
    type: boolean?
  write_protein_description:
    doc: If set, the protein description is stored as well.
    type: boolean?
  keep_unreferenced_proteins:
    doc: If set, protein hits which are not referenced by any peptide are kept.
    type: boolean?
  unmatched_action:
    doc: "If peptide sequences cannot be matched to any protein: 1) raise an error; 2) warn (unmatched PepHits will miss target/decoy annotation with downstream problems); 3) remove the hit."
    type: string?
  aaa_max:
    doc: Maximal number of ambiguous amino acids (AAAs) allowed when matching to a protein database with AAAs. AAAs are 'B', 'J', 'Z' and 'X'.
    type: long?
  mismatches_max:
    doc: Maximal number of mismatched (mm) amino acids allowed when matching to a protein database. The required runtime is exponential in the number of mm's; apply with care. MM's are allowed in addition to AAA's.
    type: long?
  IL_equivalent:
    doc: Treat the isobaric amino acids isoleucine ('I') and leucine ('L') as equivalent (indistinguishable). Also occurrences of 'J' will be treated as 'I' thus avoiding ambiguous matching.
    type: boolean?
  allow_nterm_protein_cleavage:
    doc: Allow the protein N-terminus amino acid to clip.
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
  enzyme__name:
    doc: "Enzyme which determines valid cleavage sites - e.g. trypsin cleaves after lysine (K) or arginine (R), but not before proline (P). Default: deduce from input"
    type: string?
  enzyme__specificity:
    doc: "Specificity of the enzyme. Default: deduce from input.\n  'full': both internal cleavage sites must match.\n  'semi': one of two internal cleavage sites must match.\n  'none': allow all peptide hits no matter their context (enzyme is irrelevant)."
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: PeptideIndexer
doc: Refreshes the protein references for all peptide hits.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - PeptideIndexer
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json