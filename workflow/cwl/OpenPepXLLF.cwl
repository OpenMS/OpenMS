inputs:
  in:
    doc: Input file containing the spectra.
    type: File
  database:
    doc: Input file containing the protein database.
    type: File
  decoy_database:
    doc: Input file containing the decoy protein database. Decoys can also be included in the normal database file instead (or additionally).
    type: File?
  decoy_string:
    doc: String that was appended (or prefixed - see 'prefix' flag below) to the accessions in the protein database to indicate decoy proteins.
    type: string?
  decoy_prefix:
    doc: Set to true, if the decoy_string is a prefix of accessions in the protein database. Otherwise it is a suffix.
    type: string?
  out_idXML:
    doc: Results in idXML format (at least one of these output parameters should be set, otherwise you will not have any results).
    type: string
  out_mzIdentML:
    doc: Results in mzIdentML (.mzid) format (at least one of these output parameters should be set, otherwise you will not have any results)
    type: string
  out_xquestxml:
    doc: Results in the xquest.xml format (at least one of these output parameters should be set, otherwise you will not have any results).
    type: string
  out_xquest_specxml:
    doc: Matched spectra in the xQuest .spec.xml format for spectra visualization in the xQuest results manager.
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
    doc: Width of precursor mass tolerance window
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
  precursor__corrections:
    doc: Monoisotopic peak correction. Matches candidates for possible monoisotopic precursor peaks for experimental mass m and given numbers n at masses (m - n * (C13-C12)). These should be ordered from more extreme to less extreme corrections. Numbers later in the list will be preferred in case of ambiguities.
    type: long[]?
  fragment__mass_tolerance:
    doc: Fragment mass tolerance
    type: double?
  fragment__mass_tolerance_xlinks:
    doc: Fragment mass tolerance for cross-link ions
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
  peptide__missed_cleavages:
    doc: Number of missed cleavages.
    type: long?
  peptide__enzyme:
    doc: The enzyme used for peptide digestion.
    type: string?
  cross_linker__residue1:
    doc: Comma separated residues, that the first side of a bifunctional cross-linker can attach to
    type: string[]?
  cross_linker__residue2:
    doc: Comma separated residues, that the second side of a bifunctional cross-linker can attach to
    type: string[]?
  cross_linker__mass:
    doc: Mass of the light cross-linker, linking two residues on one or two peptides
    type: double?
  cross_linker__mass_mono_link:
    doc: Possible masses of the linker, when attached to only one peptide
    type: double[]?
  cross_linker__name:
    doc: Name of the searched cross-link, used to resolve ambiguity of equal masses (e.g. DSS or BS3)
    type: string?
  algorithm__number_top_hits:
    doc: Number of top hits reported for each spectrum pair
    type: long?
  algorithm__deisotope:
    doc: Set to true, if the input spectra should be deisotoped before any other processing steps. If set to auto the spectra will be deisotoped, if the fragment mass tolerance is < 0.1 Da or < 100 ppm (0.1 Da at a mass of 1000)
    type: string?
  algorithm__use_sequence_tags:
    doc: Use sequence tags (de novo sequencing of short fragments) to filter out candidates before scoring. This will make the search faster, but can impact the sensitivity positively or negatively, depending on the dataset.
    type: boolean?
  algorithm__sequence_tag_min_length:
    doc: Minimal length of sequence tags to use for filtering candidates. Longer tags will make the search faster but much less sensitive. Ignored if 'algorithm:use_sequence_tags' is false.
    type: long?
  ions__b_ions:
    doc: Search for peaks of b-ions.
    type: string?
  ions__y_ions:
    doc: Search for peaks of y-ions.
    type: string?
  ions__a_ions:
    doc: Search for peaks of a-ions.
    type: boolean?
  ions__x_ions:
    doc: Search for peaks of x-ions.
    type: boolean?
  ions__c_ions:
    doc: Search for peaks of c-ions.
    type: boolean?
  ions__z_ions:
    doc: Search for peaks of z-ions.
    type: boolean?
  ions__neutral_losses:
    doc: Search for neutral losses of H2O and H3N.
    type: string?
outputs:
  out_idXML:
    type: File
    outputBinding:
      glob: $(inputs.out_idXML)
  out_mzIdentML:
    type: File
    outputBinding:
      glob: $(inputs.out_mzIdentML)
  out_xquestxml:
    type: File
    outputBinding:
      glob: $(inputs.out_xquestxml)
  out_xquest_specxml:
    type: File
    outputBinding:
      glob: $(inputs.out_xquest_specxml)
label: OpenPepXLLF
doc: Tool for protein-protein cross linking with label-free linkers.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - OpenPepXLLF
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json