inputs:
  in:
    doc: input file
    type: File[]
  out:
    doc: output file
    type: string
  rt_delta:
    doc: "[idXML input only] Maximum allowed retention time deviation between identifications belonging to the same spectrum."
    type: double?
  mz_delta:
    doc: "[idXML input only] Maximum allowed precursor m/z deviation between identifications belonging to the same spectrum."
    type: double?
  per_spectrum:
    doc: (only idXML) if set, mapping will be done based on exact matching of originating mzml file and spectrum_ref
    type: boolean?
  algorithm:
    doc: "Algorithm used for consensus scoring.\n* PEPMatrix: Scoring based on posterior error probabilities (PEPs) and peptide sequence similarities (scored by a substitution matrix). Requires PEPs as scores.\n* PEPIons: Scoring based on posterior error probabilities (PEPs) and fragment ion similarities ('shared peak count'). Requires PEPs as scores.\n* best: For each peptide ID, use the best score of any search engine as the consensus score. Requires the same score type in all ID runs.\n* worst: For each peptide ID, use the worst score of any search engine as the consensus score. Requires the same score type in all ID runs.\n* average:  For each peptide ID, use the average score of all search engines as the consensus. Requires the same score type in all ID runs.\n* ranks: Calculates a consensus score based on the ranks of peptide IDs in the results of different search engines. The final score is in the range (0, 1], with 1 being the best score. No requirements about score types."
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
  filter__considered_hits:
    doc: The number of top hits in each ID run that are considered for consensus scoring ('0' for all hits).
    type: long?
  filter__min_support:
    doc: For each peptide hit from an ID run, the fraction of other ID runs that must support that hit (otherwise it is removed).
    type: double?
  filter__count_empty:
    doc: Count empty ID runs (i.e. those containing no peptide hit for the current spectrum) when calculating 'min_support'?
    type: boolean?
  filter__keep_old_scores:
    doc: if set, keeps the original scores as user params
    type: boolean?
  PEPIons__mass_tolerance:
    doc: Maximum difference between fragment masses (in Da) for fragments to be considered 'shared' between peptides .
    type: double?
  PEPIons__min_shared:
    doc: The minimal number of 'shared' fragments (between two suggested peptides) that is necessary to evaluate the similarity based on shared peak count (SPC).
    type: long?
  PEPMatrix__matrix:
    doc: Substitution matrix to use for alignment-based similarity scoring
    type: string?
  PEPMatrix__penalty:
    doc: Alignment gap penalty (the same value is used for gap opening and extension)
    type: long?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: ConsensusID
doc: Computes a consensus of peptide identifications of several identification engines.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - ConsensusID
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json