inputs:
  in:
    doc: "Input file: spectra"
    type: File
  database:
    doc: "Input file: sequence database. Required unless 'digest' is set."
    type: File?
  digest:
    doc: "Input file: pre-digested sequence database. Can be used instead of 'database'. Sets all 'oligo:...' parameters."
    type: File?
  out:
    doc: "Output file: mzTab"
    type: string
  id_out:
    doc: "Output file: idXML (for visualization in TOPPView)"
    type: string
  db_out:
    doc: "Output file: oms (SQLite database)"
    type: string
  digest_out:
    doc: "Output file: sequence database digest. Ignored if 'digest' input is used."
    type: string
  lfq_out:
    doc: "Output file: targets for label-free quantification using FeatureFinderMetaboIdent ('id' input)"
    type: string
  theo_ms2_out:
    doc: "Output file: theoretical MS2 spectra for precursor mass matches"
    type: string
  exp_ms2_out:
    doc: "Output file: experimental MS2 spectra for precursor mass matches"
    type: string
  decharge_ms2:
    doc: Decharge the MS2 spectra for scoring
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
  precursor__mass_tolerance:
    doc: Precursor mass tolerance (+/- around uncharged precursor mass)
    type: double?
  precursor__mass_tolerance_unit:
    doc: Unit of precursor mass tolerance
    type: string?
  precursor__min_charge:
    doc: Minimum precursor charge to be considered
    type: long?
  precursor__max_charge:
    doc: Maximum precursor charge to be considered
    type: long?
  precursor__include_unknown_charge:
    doc: Include MS2 spectra with unknown precursor charge - try to match them in any possible charge between 'min_charge' and 'max_charge', at the risk of a higher error rate
    type: boolean?
  precursor__use_avg_mass:
    doc: Use average instead of monoisotopic precursor masses (appropriate for low-resolution instruments)
    type: boolean?
  precursor__use_adducts:
    doc: Consider possible salt adducts (see 'precursor:potential_adducts') when matching precursor masses
    type: boolean?
  precursor__potential_adducts:
    doc: "Adducts considered to explain mass differences. Format: 'Element:Charge(+/-)', i.e. the number of '+' or '-' indicates the charge, e.g. 'Ca:++' indicates +2. Only used if 'precursor:use_adducts' is set."
    type: string[]?
  precursor__isotopes:
    doc: "Correct for mono-isotopic peak misassignments. E.g.: 1 = precursor may be misassigned to the first isotopic peak. Ignored if 'use_avg_mass' is set."
    type: long[]?
  fragment__mass_tolerance:
    doc: Fragment mass tolerance (+/- around fragment m/z)
    type: double?
  fragment__mass_tolerance_unit:
    doc: Unit of fragment mass tolerance
    type: string?
  fragment__ions:
    doc: Fragment ions to include in theoretical spectra
    type: string[]?
  modifications__variable:
    doc: Variable modifications
    type: string[]?
  modifications__variable_max_per_oligo:
    doc: Maximum number of residues carrying a variable modification per candidate oligonucleotide
    type: long?
  modifications__resolve_ambiguities:
    doc: "Attempt to resolve ambiguous modifications (e.g. 'mA?' for 'mA'/'Am') based on a-B ions.\nThis incurs a performance cost because two modifications have to be considered for each case.\nRequires a-B ions to be enabled in parameter 'fragment:ions'."
    type: boolean?
  oligo__min_size:
    doc: Minimum size an oligonucleotide must have after digestion to be considered in the search
    type: long?
  oligo__max_size:
    doc: Maximum size an oligonucleotide must have after digestion to be considered in the search, leave at 0 for no limit
    type: long?
  oligo__missed_cleavages:
    doc: Number of missed cleavages
    type: long?
  oligo__enzyme:
    doc: The enzyme used for RNA digestion
    type: string?
  report__top_hits:
    doc: Maximum number of top-scoring hits per spectrum that are reported ('0' for all hits)
    type: long?
  fdr__decoy_pattern:
    doc: String used as part of the accession to annotate decoy sequences (e.g. 'DECOY_'). Leave empty to skip the FDR/q-value calculation.
    type: string?
  fdr__cutoff:
    doc: Cut-off for FDR filtering; search hits with higher q-values will be removed
    type: double?
  fdr__remove_decoys:
    doc: Do not score hits to decoy sequences and remove them when filtering
    type: boolean?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  id_out:
    type: File
    outputBinding:
      glob: $(inputs.id_out)
  db_out:
    type: File
    outputBinding:
      glob: $(inputs.db_out)
  digest_out:
    type: File
    outputBinding:
      glob: $(inputs.digest_out)
  lfq_out:
    type: File
    outputBinding:
      glob: $(inputs.lfq_out)
  theo_ms2_out:
    type: File
    outputBinding:
      glob: $(inputs.theo_ms2_out)
  exp_ms2_out:
    type: File
    outputBinding:
      glob: $(inputs.exp_ms2_out)
label: NucleicAcidSearchEngine
doc: Annotate nucleic acid identifications to MS/MS spectra.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - NucleicAcidSearchEngine
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json