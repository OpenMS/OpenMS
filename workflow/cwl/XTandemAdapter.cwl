inputs:
  in:
    doc: Input file containing MS2 spectra
    type: File
  out:
    doc: Output file containing search results
    type: string
  xml_out:
    doc: Raw output file directly from X! Tandem. Either 'out' or 'xml_out' are required. They can be used together.
    type: string
  database:
    doc: FASTA file or pro file. Non-existing relative file-names are looked up via'OpenMS.ini:id_db_dir'
    type: File
  xtandem_executable:
    doc: X! Tandem executable. Provide a full or relative path, or make sure it can be found in your PATH environment.
    type: File
  default_config_file:
    doc: Default X! Tandem configuration file. All parameters of this adapter take precedence over the file - use it for parameters not available here. A template file can be found at 'OpenMS/share/CHEMISTRY/XTandem_default_config.xml'.
    type: File?
  ignore_adapter_param:
    doc: Set this to use the configuration given in 'default_config_file' exclusively, ignoring other parameters (apart from 'in', 'out', 'database', 'xtandem_executable') set via this adapter.
    type: boolean?
  precursor_mass_tolerance:
    doc: Precursor mass tolerance
    type: double?
  fragment_mass_tolerance:
    doc: Fragment mass error
    type: double?
  precursor_error_units:
    doc: Parent monoisotopic mass error units
    type: string?
  fragment_error_units:
    doc: Fragment monoisotopic mass error units
    type: string?
  max_precursor_charge:
    doc: Maximum precursor charge ('0' to use X! Tandem default)
    type: long?
  no_isotope_error:
    doc: By default, misassignment to the first and second isotopic 13C peak are also considered. Set this flag to disable.
    type: boolean?
  fixed_modifications:
    doc: Fixed modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
    type: string[]?
  variable_modifications:
    doc: Variable modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
    type: string[]?
  minimum_fragment_mz:
    doc: Minimum fragment m/z
    type: double?
  enzyme:
    doc: The enzyme used for peptide digestion.
    type: string?
  missed_cleavages:
    doc: Number of possible cleavage sites missed by the enzyme
    type: long?
  semi_cleavage:
    doc: Require only peptide end to have a valid cleavage site, not both.
    type: boolean?
  output_results:
    doc: Which hits should be reported. All, valid ones (passing the E-Value threshold), or stochastic (failing the threshold)
    type: string?
  max_valid_expect:
    doc: Maximal E-Value of a hit to be reported (only evaluated if 'output_result' is 'valid' or 'stochastic')
    type: double?
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
  xml_out:
    type: File
    outputBinding:
      glob: $(inputs.xml_out)
label: XTandemAdapter
doc: Annotates MS/MS spectra using X! Tandem.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - XTandemAdapter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json