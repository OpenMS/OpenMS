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
  Search__enzyme:
    doc: The enzyme used for peptide digestion.
    type: string?
  Search__decoys:
    doc: Should decoys be generated?
    type: boolean?
  Search__precursor__mass_tolerance:
    doc: +/- tolerance for precursor mass.
    type: double?
  Search__precursor__mass_tolerance_unit:
    doc: Unit of precursor mass tolerance.
    type: string?
  Search__precursor__min_charge:
    doc: Minimum precursor charge to be considered.
    type: long?
  Search__precursor__max_charge:
    doc: Maximum precursor charge to be considered.
    type: long?
  Search__precursor__isotopes:
    doc: "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak)"
    type: long[]?
  Search__fragment__mass_tolerance:
    doc: Fragment mass tolerance
    type: double?
  Search__fragment__mass_tolerance_unit:
    doc: Unit of fragment m
    type: string?
  Search__modifications__fixed:
    doc: Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'
    type: string[]?
  Search__modifications__variable:
    doc: Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'
    type: string[]?
  Search__modifications__variable_max_per_peptide:
    doc: Maximum number of residues carrying a variable modification per candidate peptide
    type: long?
  Search__annotate__PSM:
    doc: Annotations added to each PSM.
    type: string[]?
  Search__peptide__min_size:
    doc: Minimum size a peptide must have after digestion to be considered in the search.
    type: long?
  Search__peptide__max_size:
    doc: Maximum size a peptide must have after digestion to be considered in the search (0 = disabled).
    type: long?
  Search__peptide__missed_cleavages:
    doc: Number of missed cleavages.
    type: long?
  Search__peptide__motif:
    doc: If set, only peptides that contain this motif (provided as RegEx) will be considered.
    type: string?
  Search__report__top_hits:
    doc: Maximum number of top scoring hits per spectrum that are reported.
    type: long?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: SimpleSearchEngine
doc: Annotates MS/MS spectra using SimpleSearchEngine.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - SimpleSearchEngine
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json