inputs:
  in:
    doc: "Input: identification results"
    type: File[]
  exp_design:
    doc: "(Currently unused) Input: experimental design"
    type: File?
  out:
    doc: "Output: identification results with scored/grouped proteins"
    type: string
  out_type:
    doc: "Output type: auto detected by file extension but can be overwritten here."
    type: string?
  protein_fdr:
    doc: Additionally calculate the target-decoy FDR on protein-level based on the posteriors
    type: boolean?
  conservative_fdr:
    doc: Use (D+1)/(T) instead of (D+1)/(T+D) for reporting protein FDRs.
    type: string?
  picked_fdr:
    doc: Use picked protein FDRs.
    type: string?
  picked_decoy_string:
    doc: If using picked protein FDRs, which decoy string was used? Leave blank for auto-detection.
    type: string?
  picked_decoy_prefix:
    doc: If using picked protein FDRs, was the decoy string a prefix or suffix? Ignored during auto-detection.
    type: string?
  greedy_group_resolution:
    doc: Post-process inference output with greedy resolution of shared peptides based on the parent protein probabilities. Also adds the resolved ambiguity groups to output.
    type: string?
  min_psms_extreme_probability:
    doc: Set PSMs with probability lower than this to this minimum probability.
    type: double?
  max_psms_extreme_probability:
    doc: Set PSMs with probability higher than this to this maximum probability.
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
  algorithm__psm_probability_cutoff:
    doc: Remove PSMs with probabilities less than this cutoff
    type: double?
  algorithm__top_PSMs:
    doc: Consider only top X PSMs per spectrum. 0 considers all.
    type: long?
  algorithm__keep_best_PSM_only:
    doc: Epifany uses the best PSM per peptide for inference. Discard the rest (true) or keepe.g. for quantification/reporting?
    type: string?
  algorithm__update_PSM_probabilities:
    doc: (Experimental:) Update PSM probabilities with their posteriors under consideration of the protein probabilities.
    type: string?
  algorithm__user_defined_priors:
    doc: (Experimental:) Uses the current protein scores as user-defined priors.
    type: boolean?
  algorithm__annotate_group_probabilities:
    doc: Annotates group probabilities for indistinguishable protein groups (indistinguishable by experimentally observed PSMs).
    type: string?
  algorithm__use_ids_outside_features:
    doc: (Only consensusXML) Also use IDs without associated features for inference?
    type: boolean?
  algorithm__model_parameters__prot_prior:
    doc: Protein prior probability ('gamma' parameter). Negative values enable grid search for this param.
    type: double?
  algorithm__model_parameters__pep_emission:
    doc: Peptide emission probability ('alpha' parameter). Negative values enable grid search for this param.
    type: double?
  algorithm__model_parameters__pep_spurious_emission:
    doc: Spurious peptide identification probability ('beta' parameter). Usually much smaller than emission from proteins. Negative values enable grid search for this param.
    type: double?
  algorithm__model_parameters__pep_prior:
    doc: Peptide prior probability (experimental, should be covered by combinations of the other params).
    type: double?
  algorithm__model_parameters__regularize:
    doc: Regularize the number of proteins that produce a peptide together (experimental, should be activated when using higher p-norms).
    type: boolean?
  algorithm__model_parameters__extended_model:
    doc: Uses information from different peptidoforms also across runs (automatically activated if an experimental design is given!)
    type: boolean?
  algorithm__loopy_belief_propagation__scheduling_type:
    doc: "(Not used yet) How to pick the next message: priority = based on difference to last message (higher = more important). fifo = first in first out. subtree = message passing follows a random spanning tree in each iteration"
    type: string?
  algorithm__loopy_belief_propagation__convergence_threshold:
    doc: Initial threshold under which MSE difference a message is considered to be converged.
    type: double?
  algorithm__loopy_belief_propagation__dampening_lambda:
    doc: Initial value for how strongly should messages be updated in each step. 0 = new message overwrites old completely (no dampening; only recommended for trees),0.5 = equal contribution of old and new message (stay below that),In-between it will be a convex combination of both. Prevents oscillations but hinders convergence.
    type: double?
  algorithm__loopy_belief_propagation__max_nr_iterations:
    doc: (Usually auto-determined by estimated but you can set a hard limit here). If not all messages converge, how many iterations should be done at max per connected component?
    type: long?
  algorithm__loopy_belief_propagation__p_norm_inference:
    doc: P-norm used for marginalization of multidimensional factors. 1 == sum-product inference (all configurations vote equally) (default),<= 0 == infinity = max-product inference (only best configurations propagate)The higher the value the more important high probability configurations get.
    type: double?
  algorithm__param_optimize__aucweight:
    doc: How important is target decoy AUC vs calibration of the posteriors? 0 = maximize calibration only, 1 = maximize AUC only, between = convex combination.
    type: double?
  algorithm__param_optimize__conservative_fdr:
    doc: Use (D+1)/(T) instead of (D+1)/(T+D) for parameter estimation.
    type: string?
  algorithm__param_optimize__regularized_fdr:
    doc: Use a regularized FDR for proteins without unique peptides.
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: Epifany
doc: Runs a Bayesian protein inference.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - Epifany
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json