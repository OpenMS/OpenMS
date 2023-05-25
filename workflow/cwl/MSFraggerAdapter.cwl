inputs:
  license:
    doc: Set to yes, if you have read and agreed to the MSFragger license terms.
    type: string
  java_executable:
    doc: The Java executable. Usually Java is on the system PATH. If Java is not found, use this parameter to specify the full path to Java
    type: File?
  java_heapmemory:
    doc: Maximum Java heap size (in MB)
    type: long?
  executable:
    doc: Path to the MSFragger executable to use; may be empty if the executable is globally available.
    type: File?
  in:
    doc: Input File with specta for MSFragger
    type: File
  out:
    doc: MSFragger output file
    type: string
  opt_out:
    doc: MSFragger optional output file
    type: string
  database:
    doc: Protein FASTA database file path
    type: File
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
  tolerance__precursor_mass_tolerance_lower:
    doc: Lower precursor mass tolerance
    type: double?
  tolerance__precursor_mass_tolerance_upper:
    doc: Upper precursor mass tolerance
    type: double?
  tolerance__precursor_mass_unit:
    doc: Unit of precursor mass tolerance
    type: string?
  tolerance__precursor_true_tolerance:
    doc: True precursor mass tolerance (window is +/- this value). Used for tie breaker of results (in spectrally ambiguous cases) and zero bin boosting in open searches (0 disables these features). This option is STRONGLY recommended for open searches.
    type: double?
  tolerance__precursor_true_unit:
    doc: Unit of precursor true tolerance
    type: string?
  tolerance__fragment_mass_tolerance:
    doc: Fragment mass tolerance (window is +/- this value)
    type: double?
  tolerance__fragment_mass_unit:
    doc: Unit of fragment mass tolerance
    type: string?
  tolerance__isotope_error:
    doc: Isotope correction for MS/MS events triggered on isotopic peaks. Should be set to 0 (disabled) for open search or 0/1/2 for correction of narrow window searches. Shifts the precursor mass window to multiples of this value multiplied by the mass of C13-C12.
    type: string?
  digest__search_enzyme_name:
    doc: Name of the enzyme to be written to the pepXML file
    type: string?
  digest__search_enzyme_cutafter:
    doc: Residues after which the enzyme cuts (specified as a string of amino acids)
    type: string?
  digest__search_enzyme_nocutbefore:
    doc: Residues that the enzyme will not cut before
    type: string?
  digest__num_enzyme_termini:
    doc: Number of enzyme termini (non-enzymatic (0), semi (1), fully (2)
    type: string?
  digest__allowed_missed_cleavage:
    doc: Allowed number of missed cleavages
    type: string?
  digest__min_length:
    doc: Minimum length of peptides to be generated during in-silico digestion
    type: long?
  digest__max_length:
    doc: Maximum length of peptides to be generated during in-silico digestion
    type: long?
  digest__mass_range_min:
    doc: Min mass of peptides to be generated (Da)
    type: double?
  digest__mass_range_max:
    doc: Max mass of peptides to be generated (Da)
    type: double?
  varmod__clip_nterm_m:
    doc: Specifies the trimming of a protein N-terminal methionine as a variable modification
    type: boolean?
  varmod__masses:
    doc: Masses for variable modifications
    type: double[]?
  varmod__syntaxes:
    doc: Syntax Strings for variable modifications
    type: string[]?
  varmod__unimod:
    doc: Variable modifications in unimod syntax, is added to mass+syntax varmod list
    type: string[]?
  varmod__enable_common:
    doc: Enable common variable modifications (15.9949 M and 42.0106 [^)
    type: boolean?
  varmod__not_allow_multiple_variable_mods_on_residue:
    doc: Do not allow any one amino acid to be modified by multiple variable modifications
    type: boolean?
  varmod__max_variable_mods_per_peptide:
    doc: Maximum total number of variable modifications per peptide
    type: string?
  varmod__max_variable_mods_combinations:
    doc: Maximum allowed number of modified variably modified peptides from each peptide sequence, (maximum of 65534). If a greater number than the maximum is generated, only the unmodified peptide is considered
    type: long?
  spectrum__minimum_peaks:
    doc: Minimum number of peaks in experimental spectrum for matching
    type: long?
  spectrum__use_topn_peaks:
    doc: Pre-process experimental spectrum to only use top N peaks
    type: long?
  spectrum__minimum_ratio:
    doc: Filters out all peaks in experimental spectrum less intense than this multiple of the base peak intensity
    type: double?
  spectrum__clear_mz_range_min:
    doc: Removes peaks in this m/z range prior to matching (minimum value). Useful for iTRAQ/TMT experiments (i.e. 0.0 150.0)
    type: double?
  spectrum__clear_mz_range_max:
    doc: Removes peaks in this m/z range prior to matching (maximum value). Useful for iTRAQ/TMT experiments (i.e. 0.0 150.0)
    type: double?
  spectrum__max_fragment_charge:
    doc: Maximum charge state for theoretical fragments to match
    type: string?
  spectrum__override_charge:
    doc: "Ignores precursor charge and uses charge state specified in precursor_charge range (parameters: spectrum:precursor_charge_min and spectrum:precursor_charge_max)"
    type: boolean?
  spectrum__precursor_charge_min:
    doc: Min charge of precursor charge range to consider. If specified, also spectrum:override_charge must be set)
    type: long?
  spectrum__precursor_charge_max:
    doc: Max charge of precursor charge range to consider. If specified, also spectrum:override_charge must be set)
    type: long?
  search__track_zero_topn:
    doc: Track top N unmodified peptide results separately from main results internally for boosting features. Should be set to a number greater than search:output_report_topN if zero bin boosting is desired
    type: long?
  search__zero_bin_accept_expect:
    doc: Ranks a zero-bin hit above all non-zero-bin hit if it has expectation less than this value
    type: double?
  search__zero_bin_mult_expect:
    doc: Multiplies expect value of PSMs in the zero-bin during results ordering (set to less than 1 for boosting)
    type: double?
  search__add_topn_complementary:
    doc: Inserts complementary ions corresponding to the top N most intense fragments in each experimental spectrum. Useful for recovery of modified peptides near C-terminus in open search. 0 disables this option
    type: long?
  search__min_fragments_modeling:
    doc: Minimum number of matched peaks in PSM for inclusion in statistical modeling
    type: long?
  search__min_matched_fragments:
    doc: Minimum number of matched peaks for PSM to be reported. MSFragger recommends a minimum of 4 for narrow window searching and 6 for open searches
    type: long?
  search__output_report_topn:
    doc: Reports top N PSMs per input spectrum
    type: long?
  search__output_max_expect:
    doc: Suppresses reporting of PSM if top hit has expectation greater than this threshold
    type: double?
  search__localize_delta_mass:
    doc: Include fragment ions mass-shifted by unknown modifications (recommended for open and mass offset searches) (0 for OFF, 1 for ON)
    type: long?
  statmod__add_cterm_peptide:
    doc: Statically add mass in Da to C-terminal of peptide
    type: double?
  statmod__add_nterm_peptide:
    doc: Statically add mass in Da to N-terminal of peptide
    type: double?
  statmod__add_cterm_protein:
    doc: Statically add mass in Da to C-terminal of protein
    type: double?
  statmod__add_nterm_protein:
    doc: Statically add mass in Da to N-terminal of protein
    type: double?
  statmod__add_G_glycine:
    doc: Statically add mass to glycine
    type: double?
  statmod__add_A_alanine:
    doc: Statically add mass to alanine
    type: double?
  statmod__add_S_serine:
    doc: Statically add mass to serine
    type: double?
  statmod__add_P_proline:
    doc: Statically add mass to proline
    type: double?
  statmod__add_V_valine:
    doc: Statically add mass to valine
    type: double?
  statmod__add_T_threonine:
    doc: Statically add mass to threonine
    type: double?
  statmod__add_C_cysteine:
    doc: Statically add mass to cysteine
    type: double?
  statmod__add_L_leucine:
    doc: Statically add mass to leucine
    type: double?
  statmod__add_I_isoleucine:
    doc: Statically add mass to isoleucine
    type: double?
  statmod__add_N_asparagine:
    doc: Statically add mass to asparagine
    type: double?
  statmod__add_D_aspartic_acid:
    doc: Statically add mass to aspartic_acid
    type: double?
  statmod__add_Q_glutamine:
    doc: Statically add mass to glutamine
    type: double?
  statmod__add_K_lysine:
    doc: Statically add mass to lysine
    type: double?
  statmod__add_E_glutamic_acid:
    doc: Statically add mass to glutamic_acid
    type: double?
  statmod__add_M_methionine:
    doc: Statically add mass to methionine
    type: double?
  statmod__add_H_histidine:
    doc: Statically add mass to histidine
    type: double?
  statmod__add_F_phenylalanine:
    doc: Statically add mass to phenylalanine
    type: double?
  statmod__add_R_arginine:
    doc: Statically add mass to arginine
    type: double?
  statmod__add_Y_tyrosine:
    doc: Statically add mass to tyrosine
    type: double?
  statmod__add_W_tryptophan:
    doc: Statically add mass to tryptophan
    type: double?
  statmod__unimod:
    doc: Fixed modifications in unimod syntax if specific mass is unknown, e.g. Carbamidomethylation (C). When multiple different masses are given for one aminoacid this parameter (unimod) will have priority.
    type: string[]?
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
  opt_out:
    type: File
    outputBinding:
      glob: $(inputs.opt_out)
label: MSFraggerAdapter
doc: "Peptide Identification with MSFragger.\nImportant note:\nThe Regents of the University of Michigan (“Michigan”) grants us permission to redistribute    \nthe MS Fragger application developed by Michigan within the OpenMS Pipeline and make available \nfor use on related service offerings supported by the University of Tubingen and the Center for\nIntegrative Bioinformatics.                                                                    \nPer the license agreement the use of the pipeline and associated materials is for academic     \nresearch, non-commercial or educational purposes. Any commercial use inquiries                 \nmust be directed to the University of Michigan Technology Transfer Office at                   \ntechtransfer@umich.edu. All right title and interest in MS Fragger shall remain with the       \nUniversity of Michigan.\n\nFor details, please see the supplied license file or                                           \nhttps://raw.githubusercontent.com/OpenMS/THIRDPARTY/master/All/MSFragger/License.txt           \n"
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MSFraggerAdapter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json