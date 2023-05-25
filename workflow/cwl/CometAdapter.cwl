inputs:
  in:
    doc: Input file
    type: File
  out:
    doc: Output file
    type: string
  database:
    doc: FASTA file
    type: File
  comet_executable:
    doc: The Comet executable. Provide a full or relative path, or make sure it can be found in your PATH environment.
    type: File
  pin_out:
    doc: Output file - for Percolator input
    type: string
  default_params_file:
    doc: Default Comet params file. All parameters of this take precedence. A template file can be generated using 'comet.exe -p'
    type: File?
  precursor_mass_tolerance:
    doc: "Precursor monoisotopic mass tolerance (Comet parameter: peptide_mass_tolerance).  See also precursor_error_units to set the unit."
    type: double?
  precursor_error_units:
    doc: "Unit of precursor monoisotopic mass tolerance for parameter precursor_mass_tolerance (Comet parameter: peptide_mass_units)"
    type: string?
  isotope_error:
    doc: This parameter controls whether the peptide_mass_tolerance takes into account possible isotope errors in the precursor mass measurement. Use -8/-4/0/4/8 only for SILAC.
    type: string?
  fragment_mass_tolerance:
    doc: "This is half the bin size, which is used to segment the MS/MS spectrum. Thus, the value should be a bit higher than for other search engines, since the bin might not be centered around the peak apex (see 'fragment_bin_offset').CAUTION: Low tolerances have heavy impact on RAM usage (since Comet uses a lot of bins in this case). Consider using use_sparse_matrix and/or spectrum_batch_size."
    type: double?
  fragment_error_units:
    doc: Fragment monoisotopic mass error units
    type: string?
  fragment_bin_offset:
    doc: "Offset of fragment bins. Recommended by Comet: low-res: 0.4, high-res: 0.0"
    type: double?
  instrument:
    doc: "Comets theoretical_fragment_ions parameter: theoretical fragment ion peak representation, high-res: sum of intensities plus flanking bins, ion trap (low-res) ms/ms: sum of intensities of central M bin only"
    type: string?
  use_A_ions:
    doc: use A ions for PSM
    type: boolean?
  use_B_ions:
    doc: use B ions for PSM
    type: string?
  use_C_ions:
    doc: use C ions for PSM
    type: boolean?
  use_X_ions:
    doc: use X ions for PSM
    type: boolean?
  use_Y_ions:
    doc: use Y ions for PSM
    type: string?
  use_Z_ions:
    doc: use Z ions for PSM
    type: boolean?
  use_NL_ions:
    doc: use neutral loss (NH3, H2O) ions from b/y for PSM
    type: boolean?
  enzyme:
    doc: The enzyme used for peptide digestion.
    type: string?
  second_enzyme:
    doc: Additional enzyme used for peptide digestion.
    type: string?
  num_enzyme_termini:
    doc: Specify the termini where the cleavage rule has to match
    type: string?
  missed_cleavages:
    doc: Number of possible cleavage sites missed by the enzyme. It has no effect if enzyme is unspecific cleavage.
    type: long?
  min_peptide_length:
    doc: Minimum peptide length to consider.
    type: long?
  max_peptide_length:
    doc: Maximum peptide length to consider.
    type: long?
  num_hits:
    doc: Number of peptide hits in output file
    type: long?
  precursor_charge:
    doc: "Precursor charge range to search (if spectrum is not annotated with a charge or if override_charge!=keep any known): 0:[num] == search all charges, 2:6 == from +2 to +6, 3:3 == +3"
    type: string?
  override_charge:
    doc: "_keep any known_: keep any precursor charge state (from input), _ignore known_: ignore known precursor charge state and use precursor_charge parameter, _ignore outside range_: ignore precursor charges outside precursor_charge range, _keep known search unknown_: keep any known precursor charge state. For unknown charge states, search as singly charged if there is no signal above the precursor m/z or use the precursor_charge range"
    type: string?
  ms_level:
    doc: MS level to analyze, valid are levels 2 (default) or 3
    type: long?
  activation_method:
    doc: If not ALL, only searches spectra of the given method
    type: string?
  digest_mass_range:
    doc: MH+ peptide mass range to analyze
    type: string?
  max_fragment_charge:
    doc: Set maximum fragment charge state to analyze as long as still lower than precursor charge - 1. (Allowed max 5)
    type: long?
  max_precursor_charge:
    doc: set maximum precursor charge state to analyze (allowed max 9)
    type: long?
  clip_nterm_methionine:
    doc: If set to true, also considers the peptide sequence w/o N-term methionine separately and applies appropriate N-term mods to it
    type: boolean?
  spectrum_batch_size:
    doc: max. number of spectra to search at a time; use 0 to search the entire scan range in one batch
    type: long?
  mass_offsets:
    doc: One or more mass offsets to search (values subtracted from deconvoluted precursor mass). Has to include 0.0 if you want the default mass to be searched.
    type: double[]?
  minimum_peaks:
    doc: Required minimum number of peaks in spectrum to search (default 10)
    type: long?
  minimum_intensity:
    doc: Minimum intensity value to read in
    type: double?
  remove_precursor_peak:
    doc: no = no removal, yes = remove all peaks around precursor m/z, charge_reduced = remove all charge reduced precursor peaks (for ETD/ECD). phosphate_loss = remove the HPO3 (-80) and H3PO4 (-98) precursor phosphate neutral loss peaks. See also remove_precursor_tolerance
    type: string?
  remove_precursor_tolerance:
    doc: one-sided tolerance for precursor removal in Thompson
    type: double?
  clear_mz_range:
    doc: for iTRAQ/TMT type data; will clear out all peaks in the specified m/z range, if not 0:0
    type: string?
  fixed_modifications:
    doc: Fixed modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
    type: string[]?
  variable_modifications:
    doc: Variable modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
    type: string[]?
  binary_modifications:
    doc: "List of modification group indices. Indices correspond to the binary modification index used by comet to group individually searched lists of variable modifications.\nNote: if set, both variable_modifications and binary_modifications need to have the same number of entries as the N-th entry corresponds to the N-th variable_modification.\n      if left empty (default), all entries are internally set to 0 generating all permutations of modified and unmodified residues.\n      For a detailed explanation please see the parameter description in the Comet help."
    type: long[]?
  max_variable_mods_in_peptide:
    doc: Set a maximum number of variable modifications per peptide
    type: long?
  require_variable_mod:
    doc: If true, requires at least one variable modification per peptide
    type: boolean?
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
  pin_out:
    type: File
    outputBinding:
      glob: $(inputs.pin_out)
label: CometAdapter
doc: Annotates MS/MS spectra using Comet.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - CometAdapter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json