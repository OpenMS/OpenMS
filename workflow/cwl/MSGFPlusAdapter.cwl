inputs:
  in:
    doc: Input file (MS-GF+ parameter '-s')
    type: File
  out:
    doc: Output file
    type: string
  mzid_out:
    doc: "Alternative output file (MS-GF+ parameter '-o')\nEither 'out' or 'mzid_out' are required. They can be used together."
    type: string
  executable:
    doc: The MSGFPlus Java archive file. Provide a full or relative path, or make sure it can be found in your PATH environment.
    type: File
  database:
    doc: Protein sequence database (FASTA file; MS-GF+ parameter '-d'). Non-existing relative filenames are looked up via 'OpenMS.ini:id_db_dir'.
    type: File
  precursor_mass_tolerance:
    doc: Precursor monoisotopic mass tolerance (MS-GF+ parameter '-t')
    type: double?
  precursor_error_units:
    doc: Unit of precursor mass tolerance (MS-GF+ parameter '-t')
    type: string?
  isotope_error_range:
    doc: Range of allowed isotope peak errors (MS-GF+ parameter '-ti'). Takes into account the error introduced by choosing a non-monoisotopic peak for fragmentation. Combined with 'precursor_mass_tolerance'/'precursor_error_units', this determines the actual precursor mass tolerance. E.g. for experimental mass 'exp' and calculated mass 'calc', '-precursor_mass_tolerance 20 -precursor_error_units ppm -isotope_error_range -1,2' tests '|exp - calc - n * 1.00335 Da| < 20 ppm' for n = -1, 0, 1, 2.
    type: string?
  fragment_method:
    doc: Fragmentation method ('from_spectrum' relies on spectrum meta data and uses CID as fallback option; MS-GF+ parameter '-m')
    type: string?
  instrument:
    doc: Instrument that generated the data ('low_res'/'high_res' refer to LCQ and LTQ instruments; MS-GF+ parameter '-inst')
    type: string?
  enzyme:
    doc: "Enzyme used for digestion, or type of cleavage. Note: MS-GF+ does not support blocking rules. (MS-GF+ parameter '-e')"
    type: string?
  protocol:
    doc: Labeling or enrichment protocol used, if any (MS-GF+ parameter '-p')
    type: string?
  tryptic:
    doc: Level of cleavage specificity required (MS-GF+ parameter '-ntt')
    type: string?
  min_precursor_charge:
    doc: Minimum precursor ion charge (only used for spectra without charge information; MS-GF+ parameter '-minCharge')
    type: long?
  max_precursor_charge:
    doc: Maximum precursor ion charge (only used for spectra without charge information; MS-GF+ parameter '-maxCharge')
    type: long?
  min_peptide_length:
    doc: Minimum peptide length to consider (MS-GF+ parameter '-minLength')
    type: long?
  max_peptide_length:
    doc: Maximum peptide length to consider (MS-GF+ parameter '-maxLength')
    type: long?
  matches_per_spec:
    doc: Number of matches per spectrum to be reported (MS-GF+ parameter '-n')
    type: long?
  add_features:
    doc: Output additional features (MS-GF+ parameter '-addFeatures'). This is required by Percolator and hence by default enabled.
    type: string?
  max_mods:
    doc: Maximum number of modifications per peptide. If this value is large, the search may take very long.
    type: long?
  max_missed_cleavages:
    doc: "Maximum number of missed cleavages allowed for a peptide to be considered for scoring. (default: -1 meaning unlimited)"
    type: long?
  tasks:
    doc: "(Override the number of tasks to use on the threads; Default: (internally calculated based on inputs))\n   More tasks than threads will reduce the memory requirements of the search, but will be slower (how much depends on the inputs).\n   1 <= tasks <= numThreads: will create one task per thread, which is the original behavior.\n   tasks = 0: use default calculation - minimum of: (threads*3) and (numSpectra/250).\n   tasks < 0: multiply number of threads by abs(tasks) to determine number of tasks (i.e., -2 means \"2 * numThreads\" tasks).\n   One task per thread will use the most memory, but will usually finish the fastest.\n   2-3 tasks per thread will use comparably less memory, but may cause the search to take 1.5 to 2 times as long."
    type: long?
  fixed_modifications:
    doc: Fixed modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
    type: string[]?
  variable_modifications:
    doc: Variable modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
    type: string[]?
  legacy_conversion:
    doc: Use the indirect conversion of MS-GF+ results to idXML via export to TSV. Try this only if the default conversion takes too long or uses too much memory.
    type: boolean?
  conf:
    doc: Optional MSGF+ configuration file (passed as -conf <file> to MSGF+). See documentation for examples. Parameters of the adapter take precedence. Use conf file only for settings not available here (for example, any fixed/var modifications, in the conf file will be ignored, since they are provided via -mod flag)
    type: File?
  java_executable:
    doc: The Java executable. Usually Java is on the system PATH. If Java is not found, use this parameter to specify the full path to Java
    type: File?
  java_memory:
    doc: Maximum Java heap size (in MB)
    type: long?
  java_permgen:
    doc: Maximum Java permanent generation space (in MB); only for Java 7 and below
    type: long?
  reindex:
    doc: Recalculate peptide to protein association using OpenMS. Annotates target-decoy information.
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
  PeptideIndexing__decoy_string:
    doc: String that was appended (or prefixed - see 'decoy_string_position' flag below) to the accessions in the protein database to indicate decoy proteins. If empty (default), it's determined automatically (checking for common terms, both as prefix and suffix).
    type: string?
  PeptideIndexing__decoy_string_position:
    doc: Is the 'decoy_string' prepended (prefix) or appended (suffix) to the protein accession? (ignored if decoy_string is empty)
    type: string?
  PeptideIndexing__missing_decoy_action:
    doc: "Action to take if NO peptide was assigned to a decoy protein (which indicates wrong database or decoy string): 'error' (exit with error, no output), 'warn' (exit with success, warning message), 'silent' (no action is taken, not even a warning)"
    type: string?
  PeptideIndexing__write_protein_sequence:
    doc: If set, the protein sequences are stored as well.
    type: boolean?
  PeptideIndexing__write_protein_description:
    doc: If set, the protein description is stored as well.
    type: boolean?
  PeptideIndexing__keep_unreferenced_proteins:
    doc: If set, protein hits which are not referenced by any peptide are kept.
    type: boolean?
  PeptideIndexing__unmatched_action:
    doc: "If peptide sequences cannot be matched to any protein: 1) raise an error; 2) warn (unmatched PepHits will miss target/decoy annotation with downstream problems); 3) remove the hit."
    type: string?
  PeptideIndexing__aaa_max:
    doc: Maximal number of ambiguous amino acids (AAAs) allowed when matching to a protein database with AAAs. AAAs are 'B', 'J', 'Z' and 'X'.
    type: long?
  PeptideIndexing__mismatches_max:
    doc: Maximal number of mismatched (mm) amino acids allowed when matching to a protein database. The required runtime is exponential in the number of mm's; apply with care. MM's are allowed in addition to AAA's.
    type: long?
  PeptideIndexing__IL_equivalent:
    doc: Treat the isobaric amino acids isoleucine ('I') and leucine ('L') as equivalent (indistinguishable). Also occurrences of 'J' will be treated as 'I' thus avoiding ambiguous matching.
    type: boolean?
  PeptideIndexing__allow_nterm_protein_cleavage:
    doc: Allow the protein N-terminus amino acid to clip.
    type: string?
  PeptideIndexing__enzyme__name:
    doc: "Enzyme which determines valid cleavage sites - e.g. trypsin cleaves after lysine (K) or arginine (R), but not before proline (P). Default: deduce from input"
    type: string?
  PeptideIndexing__enzyme__specificity:
    doc: "Specificity of the enzyme. Default: deduce from input.\n  'full': both internal cleavage sites must match.\n  'semi': one of two internal cleavage sites must match.\n  'none': allow all peptide hits no matter their context (enzyme is irrelevant)."
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  mzid_out:
    type: File
    outputBinding:
      glob: $(inputs.mzid_out)
label: MSGFPlusAdapter
doc: MS/MS database search using MS-GF+.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MSGFPlusAdapter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json