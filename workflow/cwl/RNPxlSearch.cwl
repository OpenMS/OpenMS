inputs:
  in:
    doc: "input file "
    type: File
  database:
    doc: "input file "
    type: File
  out:
    doc: "output file "
    type: string
  out_tsv:
    doc: tsv output file
    type: string
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
  precursor__mass_tolerance:
    doc: Precursor mass tolerance (+/- around precursor m/z)
    type: double?
  precursor__mass_tolerance_unit:
    doc: Unit of precursor mass tolerance.
    type: string?
  precursor__min_charge:
    doc: Minimum precursor charge to be considered.
    type: long?
  precursor__max_charge:
    doc: Maximum precursor charge to be considered.
    type: long?
  precursor__isotopes:
    doc: "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak)"
    type: long[]?
  fragment__mass_tolerance:
    doc: Fragment mass tolerance (+/- around fragment m/z)
    type: double?
  fragment__mass_tolerance_unit:
    doc: Unit of fragment m
    type: string?
  modifications__fixed:
    doc: Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'
    type: string[]?
  modifications__variable:
    doc: Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'
    type: string[]?
  modifications__variable_max_per_peptide:
    doc: Maximum number of residues carrying a variable modification per candidate peptide
    type: long?
  peptide__min_size:
    doc: Minimum size a peptide must have after digestion to be considered in the search.
    type: long?
  peptide__max_size:
    doc: Maximum size a peptide may have after digestion to be considered in the search.
    type: long?
  peptide__missed_cleavages:
    doc: Number of missed cleavages.
    type: long?
  peptide__enzyme:
    doc: The enzyme used for peptide digestion.
    type: string?
  report__top_hits:
    doc: Maximum number of top scoring hits per spectrum that are reported.
    type: long?
  RNPxl__length:
    doc: Oligonucleotide maximum length. 0 = disable search for RNA variants.
    type: long?
  RNPxl__sequence:
    doc: Sequence to restrict the generation of oligonucleotide chains. (disabled for empty sequence)
    type: string?
  RNPxl__target_nucleotides:
    doc: "format:  target nucleotide=empirical formula of nucleoside monophosphate \n e.g. A=C10H14N5O7P, ..., U=C10H14N5O7P, X=C9H13N2O8PS  where X represents e.g. tU \n or e.g. Y=C10H14N5O7PS where Y represents tG"
    type: string[]?
  RNPxl__nt_groups:
    doc: "Restrict which nucleotides can cooccur in a precursor adduct to be able to search both RNA and DNA (Formate e.g.: AU CG)."
    type: string[]?
  RNPxl__mapping:
    doc: "format: source->target e.g. A->A, ..., U->U, U->X"
    type: string[]?
  RNPxl__can_cross_link:
    doc: "format: 'U' if only U forms cross-links. 'CATG' if C, A, G, and T form cross-links."
    type: string?
  RNPxl__fragment_adducts:
    doc: "format: [target nucleotide]:[formula] or [precursor adduct]->[fragment adduct formula];[name]: e.g., 'U:C9H10N2O5;U-H3PO4' or 'U:U-H2O->C9H11N2O8P1;U-H2O',"
    type: string[]?
  RNPxl__modifications:
    doc: "format: empirical formula e.g -H2O, ..., H2O+PO3"
    type: string[]?
  RNPxl__scoring:
    doc: "Scoring algorithm used in prescoring (fast: total-loss, slow: all losses)."
    type: string?
  RNPxl__decoys:
    doc: Generate decoy sequences and spectra.
    type: boolean?
  RNPxl__CysteineAdduct:
    doc: Use this flag if the +152 adduct is expected.
    type: boolean?
  RNPxl__filter_fractional_mass:
    doc: Use this flag to filter non-crosslinks by fractional mass.
    type: boolean?
  RNPxl__carbon_labeled_fragments:
    doc: Generate fragment shifts assuming full labeling of carbon (e.g. completely labeled U13).
    type: boolean?
  RNPxl__only_xl:
    doc: Only search cross-links and ignore non-cross-linked peptides.
    type: boolean?
  RNPxl__filter_small_peptide_mass:
    doc: Filter precursor that can only correspond to non-crosslinks by mass.
    type: double?
  RNPxl__marker_ions_tolerance:
    doc: Tolerance used to determine marker ions (Da).
    type: double?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_tsv:
    type: File
    outputBinding:
      glob: $(inputs.out_tsv)
label: RNPxlSearch
doc: Annotate RNA/DNA-peptide cross-links in MS/MS spectra.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - RNPxlSearch
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json