inputs:
  in:
    doc: "input file "
    type: File
  out:
    doc: "output file "
    type: string
  var_mods:
    doc: Keep only peptide hits with variable modifications (as defined in the 'SearchParameters' section of the input file).
    type: boolean?
  remove_duplicate_psm:
    doc: Removes duplicated PSMs per spectrum and retains the one with the higher score.
    type: boolean?
  remove_shared_peptides:
    doc: Only peptides matching exactly one protein are kept. Remember that isoforms count as different proteins!
    type: boolean?
  keep_unreferenced_protein_hits:
    doc: Proteins not referenced by a peptide are retained in the IDs.
    type: boolean?
  remove_decoys:
    doc: Remove proteins according to the information in the user parameters. Usually used in combination with 'delete_unreferenced_peptide_hits'.
    type: boolean?
  delete_unreferenced_peptide_hits:
    doc: Peptides not referenced by any protein are deleted in the IDs. Usually used in combination with 'score:prot' or 'thresh:prot'.
    type: boolean?
  remove_peptide_hits_by_metavalue:
    doc: Expects a 3-tuple (=3 entries in the list), i.e. <name> 'lt|eq|gt|ne' <value>; the first is the name of meta value, followed by the comparison operator (equal, less, greater, not equal) and the value to compare to. All comparisons are done after converting the given value to the corresponding data value type of the meta value (for lists, this simply compares length, not content!)!
    type: string[]?
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
  precursor__rt:
    doc: Retention time range to extract.
    type: string?
  precursor__mz:
    doc: Mass-to-charge range to extract.
    type: string?
  precursor__length:
    doc: Keep only peptide hits with a sequence length in this range.
    type: string?
  precursor__charge:
    doc: Keep only peptide hits with charge states in this range.
    type: string?
  score__pep:
    doc: The score which should be reached by a peptide hit to be kept.
    type: double?
  score__prot:
    doc: The score which should be reached by a protein hit to be kept. All proteins are filtered based on their singleton scores irrespective of grouping. Use in combination with 'delete_unreferenced_peptide_hits' to remove affected peptides.
    type: double?
  score__protgroup:
    doc: The score which should be reached by a protein group to be kept. Performs group level score filtering (including groups of single proteins). Use in combination with 'delete_unreferenced_peptide_hits' to remove affected peptides.
    type: double?
  whitelist__proteins:
    doc: "Filename of a FASTA file containing protein sequences.\nAll peptides that are not referencing a protein in this file are removed.\nAll proteins whose accessions are not present in this file are removed."
    type: File?
  whitelist__protein_accessions:
    doc: "All peptides that do not reference at least one of the provided protein accession are removed.\nOnly proteins of the provided list are retained."
    type: string[]?
  whitelist__peptides:
    doc: "Only peptides with the same sequence and modification assignment as any peptide in this file are kept. Use with 'whitelist:ignore_modifications' to only compare by sequence.\n"
    type: File?
  whitelist__ignore_modifications:
    doc: Compare whitelisted peptides by sequence only.
    type: boolean?
  whitelist__modifications:
    doc: Keep only peptides with sequences that contain (any of) the selected modification(s)
    type: string[]?
  blacklist__proteins:
    doc: "Filename of a FASTA file containing protein sequences.\nAll peptides that are referencing a protein in this file are removed.\nAll proteins whose accessions are present in this file are removed."
    type: File?
  blacklist__protein_accessions:
    doc: "All peptides that reference at least one of the provided protein accession are removed.\nOnly proteins not in the provided list are retained."
    type: string[]?
  blacklist__peptides:
    doc: "Peptides with the same sequence and modification assignment as any peptide in this file are filtered out. Use with 'blacklist:ignore_modifications' to only compare by sequence.\n"
    type: File?
  blacklist__ignore_modifications:
    doc: Compare blacklisted peptides by sequence only.
    type: boolean?
  blacklist__modifications:
    doc: Remove all peptides with sequences that contain (any of) the selected modification(s)
    type: string[]?
  blacklist__RegEx:
    doc: Remove all peptides with (unmodified) sequences matched by the RegEx e.g. [BJXZ] removes ambiguous peptides.
    type: string?
  in_silico_digestion__fasta:
    doc: fasta protein sequence database.
    type: File?
  in_silico_digestion__enzyme:
    doc: enzyme used for the digestion of the sample
    type: string?
  in_silico_digestion__specificity:
    doc: Specificity of the filter
    type: string?
  in_silico_digestion__missed_cleavages:
    doc: "range of allowed missed cleavages in the peptide sequences\nBy default missed cleavages are ignored"
    type: long?
  in_silico_digestion__methionine_cleavage:
    doc: Allow methionine cleavage at the N-terminus of the protein.
    type: boolean?
  missed_cleavages__number_of_missed_cleavages:
    doc: "range of allowed missed cleavages in the peptide sequences.\nFor example: 0:1 -> peptides with two or more missed cleavages will be removed,\n0:0 -> peptides with any missed cleavages will be removed"
    type: string?
  missed_cleavages__enzyme:
    doc: enzyme used for the digestion of the sample
    type: string?
  rt__p_value:
    doc: Retention time filtering by the p-value predicted by RTPredict.
    type: double?
  rt__p_value_1st_dim:
    doc: Retention time filtering by the p-value predicted by RTPredict for first dimension.
    type: double?
  mz__error:
    doc: Filtering by deviation to theoretical mass (disabled for negative values).
    type: double?
  mz__unit:
    doc: Absolute or relative error.
    type: string?
  best__n_spectra:
    doc: Keep only the 'n' best spectra (i.e., PeptideIdentifications) (for n > 0). A spectrum is considered better if it has a higher scoring peptide hit than the other spectrum.
    type: long?
  best__n_peptide_hits:
    doc: Keep only the 'n' highest scoring peptide hits per spectrum (for n > 0).
    type: long?
  best__n_protein_hits:
    doc: Keep only the 'n' highest scoring protein hits (for n > 0).
    type: long?
  best__strict:
    doc: "Keep only the highest scoring peptide hit.\nSimilar to n_peptide_hits=1, but if there are ties between two or more highest scoring hits, none are kept."
    type: boolean?
  best__n_to_m_peptide_hits:
    doc: Peptide hit rank range to extracts
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: IDFilter
doc: Filters results from protein or peptide identification engines based on different criteria.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - IDFilter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json