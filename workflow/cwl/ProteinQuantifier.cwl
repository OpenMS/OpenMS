inputs:
  in:
    doc: Input file
    type: File
  protein_groups:
    doc: "Protein inference results for the identification runs that were used to annotate the input (e.g. from ProteinProphet via IDFileConverter or Fido via FidoAdapter).\nInformation about indistinguishable proteins will be used for protein quantification."
    type: File?
  design:
    doc: input file containing the experimental design
    type: File?
  out:
    doc: Output file for protein abundances
    type: string
  peptide_out:
    doc: Output file for peptide abundances
    type: string
  mztab:
    doc: Output file (mzTab)
    type: string
  method:
    doc: "- top - quantify based on three most abundant peptides (number can be changed in 'top').\n- iBAQ (intensity based absolute quantification), calculate the sum of all peptide peak intensities divided by the number of theoretically observable tryptic peptides (https://rdcu.be/cND1J). Warning: only consensusXML or featureXML input is allowed!"
    type: string?
  best_charge_and_fraction:
    doc: "Distinguish between fraction and charge states of a peptide. For peptides, abundances will be reported separately for each fraction and charge;\nfor proteins, abundances will be computed based only on the most prevalent charge observed of each peptide (over all fractions).\nBy default, abundances are summed over all charge states."
    type: boolean?
  greedy_group_resolution:
    doc: Pre-process identifications with greedy resolution of shared peptides based on the protein group probabilities. (Only works with an idXML file given as protein_groups parameter).
    type: boolean?
  ratios:
    doc: "Add the log2 ratios of the abundance values to the output. Format: log_2(x_0/x_0) <sep> log_2(x_1/x_0) <sep> log_2(x_2/x_0) ..."
    type: boolean?
  ratiosSILAC:
    doc: "Add the log2 ratios for a triple SILAC experiment to the output. Only applicable to consensus maps of exactly three sub-maps. Format: log_2(heavy/light) <sep> log_2(heavy/middle) <sep> log_2(middle/light)"
    type: boolean?
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
  top__N:
    doc: Calculate protein abundance from this number of proteotypic peptides (most abundant first; '0' for all)
    type: long?
  top__aggregate:
    doc: Aggregation method used to compute protein abundances from peptide abundances
    type: string?
  top__include_all:
    doc: Include results for proteins with fewer proteotypic peptides than indicated by 'N' (no effect if 'N' is 0 or 1)
    type: boolean?
  consensus__normalize:
    doc: Scale peptide abundances so that medians of all samples are equal
    type: boolean?
  consensus__fix_peptides:
    doc: "Use the same peptides for protein quantification across all samples.\nWith 'N 0',all peptides that occur in every sample are considered.\nOtherwise ('N'), the N peptides that occur in the most samples (independently of each other) are selected,\nbreaking ties by total abundance (there is no guarantee that the best co-ocurring peptides are chosen!)."
    type: boolean?
  format__separator:
    doc: Character(s) used to separate fields; by default, the 'tab' character is used
    type: string?
  format__quoting:
    doc: "Method for quoting of strings: 'none' for no quoting, 'double' for quoting with doubling of embedded quotes,\n'escape' for quoting with backslash-escaping of embedded quotes"
    type: string?
  format__replacement:
    doc: If 'quoting' is 'none', used to replace occurrences of the separator in strings before writing
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  peptide_out:
    type: File
    outputBinding:
      glob: $(inputs.peptide_out)
  mztab:
    type: File
    outputBinding:
      glob: $(inputs.mztab)
label: ProteinQuantifier
doc: Compute peptide and protein abundances
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - ProteinQuantifier
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json