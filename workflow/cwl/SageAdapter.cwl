# Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: Apache-2.0
label: SageAdapter
doc: Annotates MS/MS spectra using Sage.
inputs:
  in:
    doc: Input files separated by blank
    type: File[]
  out:
    doc: Single output file containing all search results.
    type: string
  database:
    doc: FASTA file
    type: File
  sage_executable:
    doc: The Sage executable. Provide a full or relative path, or make sure it can be found in your PATH environment.
    type: File
  decoy_prefix:
    doc: "Prefix on protein accession used to distinguish decoy from target proteins. NOTE: Decoy suffix is currently not supported by sage."
    type: string?
  batch_size:
    doc: "Number of files to load and search in parallel (default = # of CPUs/2)"
    type: long?
  precursor_tol_left:
    doc: Start (left side) of the precursor tolerance window w.r.t. precursor location. Usually used with negative values smaller or equal to the 'right' counterpart.
    type: double?
  precursor_tol_right:
    doc: End (right side) of the precursor tolerance window w.r.t. precursor location. Usually used with positive values larger or equal to the 'left' counterpart.
    type: double?
  precursor_tol_unit:
    doc: Unit of precursor tolerance (ppm or Da)
    type: string?
  fragment_tol_left:
    doc: Start (left side) of the fragment tolerance window w.r.t. precursor location. Usually used with negative values smaller or equal to the 'right' counterpart.
    type: double?
  fragment_tol_right:
    doc: End (right side) of the fragment tolerance window w.r.t. precursor location. Usually used with positive values larger or equal to the 'left' counterpart.
    type: double?
  fragment_tol_unit:
    doc: Unit of fragment tolerance (ppm or Da)
    type: string?
  min_matched_peaks:
    doc: Minimum number of b+y ions required to match for PSM to be reported
    type: long?
  min_peaks:
    doc: Minimum number of peaks required for a spectrum to be considered
    type: long?
  max_peaks:
    doc: Take the top N most intense MS2 peaks only for matching
    type: long?
  report_psms:
    doc: Number of hits (PSMs) to report for each spectrum
    type: long?
  bucket_size:
    doc: "How many fragments are in each internal mass bucket (default: 8192 for hi-res data). Try increasing it to 32k or 64k for low-res. See also: fragment_tol_*"
    type: long?
  min_len:
    doc: Minimum peptide length
    type: long?
  max_len:
    doc: Maximum peptide length
    type: long?
  missed_cleavages:
    doc: Number of missed cleavages
    type: long?
  fragment_min_mz:
    doc: Minimum fragment m/z
    type: double?
  fragment_max_mz:
    doc: Maximum fragment m/z
    type: double?
  peptide_min_mass:
    doc: Minimum monoisotopic peptide mass to consider a peptide from the DB
    type: double?
  peptide_max_mass:
    doc: Maximum monoisotopic peptide mass to consider a peptide from the DB
    type: double?
  min_ion_index:
    doc: Minimum ion index to consider for preliminary scoring. Default = 2 to skip b1/y1 AND (sic) b2/y2 ions that are often missing.
    type: long?
  max_variable_mods:
    doc: Maximum number of variable modifications
    type: long?
  isotope_error_range:
    doc: Range of (C13) isotope errors to consider for precursor.Can be negative. E.g. '-1,3' for considering '-1/0/1/2/3'
    type: string?
  charges:
    doc: Range of precursor charges to consider if not annotated in the file.
    type: string?
  enzyme:
    doc: The enzyme used for peptide digestion.
    type: string?
  fixed_modifications:
    doc: Fixed modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
    type: string[]?
  variable_modifications:
    doc: Variable modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
    type: string[]?
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
  test:
    doc: Enables the test mode (needed for internal use only)
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
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - SageAdapter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json
