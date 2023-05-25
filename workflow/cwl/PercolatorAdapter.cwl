inputs:
  in:
    doc: Input file(s)
    type: File[]?
  in_decoy:
    doc: Input decoy file(s) in case of separate searches
    type: File[]?
  in_osw:
    doc: Input file in OSW format
    type: File?
  out:
    doc: Output file
    type: string
  out_pin:
    doc: Write pin file (e.g., for debugging)
    type: string
  out_pout_target:
    doc: Write pout file (e.g., for debugging)
    type: string
  out_pout_decoy:
    doc: Write pout file (e.g., for debugging)
    type: string
  out_pout_target_proteins:
    doc: Write pout file (e.g., for debugging)
    type: string
  out_pout_decoy_proteins:
    doc: Write pout file (e.g., for debugging)
    type: string
  out_type:
    doc: "Output file type -- default: determined from file extension or content."
    type: string?
  enzyme:
    doc: "Type of enzyme: no_enzyme,elastase,pepsin,proteinasek,thermolysin,chymotrypsin,lys-n,lys-c,arg-c,asp-n,glu-c,trypsin,trypsinp"
    type: string?
  percolator_executable:
    doc: The Percolator executable. Provide a full or relative path, or make sure it can be found in your PATH environment.
    type: File
  peptide_level_fdrs:
    doc: Calculate peptide-level FDRs instead of PSM-level FDRs.
    type: boolean?
  protein_level_fdrs:
    doc: Use the picked protein-level FDR to infer protein probabilities. Use the -fasta option and -decoy_pattern to set the Fasta file and decoy pattern.
    type: boolean?
  osw_level:
    doc: "OSW: the data level selected for scoring."
    type: string?
  score_type:
    doc: Type of the peptide main score
    type: string?
  generic_feature_set:
    doc: Use only generic (i.e. not search engine specific) features. Generating search engine specific features for common search engines by PSMFeatureExtractor will typically boost the identification rate significantly.
    type: boolean?
  subset_max_train:
    doc: Only train an SVM on a subset of <x> PSMs, and use the resulting score vector to evaluate the other PSMs. Recommended when analyzing huge numbers (>1 million) of PSMs. When set to 0, all PSMs are used for training as normal.
    type: long?
  cpos:
    doc: Cpos, penalty for mistakes made on positive examples. Set by cross validation if not specified.
    type: double?
  cneg:
    doc: Cneg, penalty for mistakes made on negative examples. Set by cross validation if not specified.
    type: double?
  testFDR:
    doc: False discovery rate threshold for evaluating best cross validation result and the reported end result.
    type: double?
  trainFDR:
    doc: False discovery rate threshold to define positive examples in training. Set to testFDR if 0.
    type: double?
  maxiter:
    doc: Maximal number of iterations
    type: long?
  nested_xval_bins:
    doc: Number of nested cross-validation bins in the 3 splits.
    type: long?
  quick_validation:
    doc: Quicker execution by reduced internal cross-validation.
    type: boolean?
  weights:
    doc: Output final weights to the given file
    type: string
  init_weights:
    doc: Read initial weights to the given file
    type: File?
  static:
    doc: Use static model (requires init-weights parameter to be set)
    type: boolean?
  default_direction:
    doc: The most informative feature given as the feature name, can be negated to indicate that a lower value is better.
    type: string?
  verbose:
    doc: "Set verbosity of output: 0=no processing info, 5=all."
    type: long?
  unitnorm:
    doc: Use unit normalization [0-1] instead of standard deviation normalization
    type: boolean?
  test_each_iteration:
    doc: Measure performance on test set each iteration
    type: boolean?
  override:
    doc: Override error check and do not fall back on default score vector in case of suspect score vector
    type: boolean?
  seed:
    doc: Setting seed of the random number generator.
    type: long?
  doc:
    doc: Include description of correct features
    type: long?
  klammer:
    doc: Retention time features calculated as in Klammer et al. Only available if -doc is set
    type: boolean?
  fasta:
    doc: Provide the fasta file as the argument to this flag, which will be used for protein grouping based on an in-silico digest (only valid if option -protein_level_fdrs is active).
    type: File?
  decoy_pattern:
    doc: Define the text pattern to identify the decoy proteins and/or PSMs, set this up if the label that identifies the decoys in the database is not the default (Only valid if option -protein_level_fdrs is active).
    type: string?
  post_processing_tdc:
    doc: Use target-decoy competition to assign q-values and PEPs.
    type: boolean?
  train_best_positive:
    doc: Enforce that, for each spectrum, at most one PSM is included in the positive set during each training iteration. If the user only provides one PSM per spectrum, this filter will have no effect.
    type: boolean?
  ipf_max_peakgroup_pep:
    doc: "OSW/IPF: Assess transitions only for candidate peak groups until maximum posterior error probability."
    type: double?
  ipf_max_transition_isotope_overlap:
    doc: "OSW/IPF: Maximum isotope overlap to consider transitions in IPF."
    type: double?
  ipf_min_transition_sn:
    doc: "OSW/IPF: Minimum log signal-to-noise level to consider transitions in IPF. Set -1 to disable this filter."
    type: double?
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
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_pin:
    type: File
    outputBinding:
      glob: $(inputs.out_pin)
  out_pout_target:
    type: File
    outputBinding:
      glob: $(inputs.out_pout_target)
  out_pout_decoy:
    type: File
    outputBinding:
      glob: $(inputs.out_pout_decoy)
  out_pout_target_proteins:
    type: File
    outputBinding:
      glob: $(inputs.out_pout_target_proteins)
  out_pout_decoy_proteins:
    type: File
    outputBinding:
      glob: $(inputs.out_pout_decoy_proteins)
  weights:
    type: File
    outputBinding:
      glob: $(inputs.weights)
label: PercolatorAdapter
doc: Facilitate input to Percolator and reintegrate.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - PercolatorAdapter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json