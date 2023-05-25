inputs:
  in:
    doc: input file in mzML format
    type: File
  out:
    doc: output file in idXML format
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
  algorithm__max_number_aa_per_decomp:
    doc: maximal amino acid frequency per decomposition
    type: long?
  algorithm__tryptic_only:
    doc: if set to true only tryptic peptides are reported
    type: string?
  algorithm__precursor_mass_tolerance:
    doc: precursor mass tolerance
    type: double?
  algorithm__fragment_mass_tolerance:
    doc: fragment mass tolerance
    type: double?
  algorithm__max_number_pivot:
    doc: maximal number of pivot ions to be used
    type: long?
  algorithm__max_subscore_number:
    doc: maximal number of solutions of a subsegment that are kept
    type: long?
  algorithm__decomp_weights_precision:
    doc: precision used to calculate the decompositions, this only affects cache usage!
    type: double?
  algorithm__double_charged_iso_threshold:
    doc: minimal isotope intensity correlation of doubly charged ions to be used to score the single scored ions
    type: double?
  algorithm__max_mz:
    doc: maximal m/z value used to calculate isotope distributions
    type: double?
  algorithm__min_mz:
    doc: minimal m/z value used to calculate the isotope distributions
    type: double?
  algorithm__max_isotope_to_score:
    doc: max isotope peak to be considered in the scoring
    type: long?
  algorithm__max_decomp_weight:
    doc: maximal m/z difference used to calculate the decompositions
    type: double?
  algorithm__max_isotope:
    doc: max isotope used in the theoretical spectra to score
    type: long?
  algorithm__missed_cleavages:
    doc: maximal number of missed cleavages allowed per peptide
    type: long?
  algorithm__number_of_hits:
    doc: maximal number of hits which are reported per spectrum
    type: long?
  algorithm__estimate_precursor_mz:
    doc: "If set to true, the precursor charge will be estimated, e.g. from the precursor peaks of the ETD spectrum.\nThe input is believed otherwise."
    type: string?
  algorithm__number_of_prescoring_hits:
    doc: how many sequences are kept after first rough scoring for better scoring
    type: long?
  algorithm__fixed_modifications:
    doc: fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
    type: string[]?
  algorithm__variable_modifications:
    doc: variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
    type: string[]?
  algorithm__residue_set:
    doc: The predefined amino acid set that should be used, see doc of ResidueDB for possible residue sets
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: CompNovo
doc: Performs a de novo peptide identification using the CompNovo engine.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - CompNovo
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json