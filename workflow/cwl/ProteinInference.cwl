inputs:
  in:
    doc: input file(s)
    type: File[]
  out:
    doc: output file
    type: string
  out_type:
    doc: output file type
    type: string?
  merge_runs:
    doc: If your idXML contains multiple runs, merge them beforehand? Otherwise performs inference separately per run.
    type: string?
  protein_fdr:
    doc: Additionally calculate the target-decoy FDR on protein-level after inference
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
  Merging__annotate_origin:
    doc: If true, adds a map_index MetaValue to the PeptideIDs to annotate the IDRun they came from.
    type: string?
  Merging__allow_disagreeing_settings:
    doc: Force merging of disagreeing runs. Use at your own risk.
    type: boolean?
  Algorithm__min_peptides_per_protein:
    doc: Minimal number of peptides needed for a protein identification. If set to zero, unmatched proteins get a score of -Infinity. If bigger than zero, proteins with less peptides are filtered and evidences removed from the PSMs. PSMs that do not reference any proteins anymore are removed but the spectrum info is kept.
    type: long?
  Algorithm__score_aggregation_method:
    doc: How to aggregate scores of peptides matching to the same protein?
    type: string?
  Algorithm__treat_charge_variants_separately:
    doc: If this is true, different charge variants of the same peptide sequence count as individual evidences.
    type: string?
  Algorithm__treat_modification_variants_separately:
    doc: If this is true, different modification variants of the same peptide sequence count as individual evidences.
    type: string?
  Algorithm__use_shared_peptides:
    doc: "If this is true, shared peptides are used as evidences. Note: shared_peptides are not deleted and potentially resolved in postprocessing as well."
    type: string?
  Algorithm__skip_count_annotation:
    doc: If this is set, peptide counts won't be annotated at the proteins.
    type: boolean?
  Algorithm__annotate_indistinguishable_groups:
    doc: If this is true, calculates and annotates indistinguishable protein groups.
    type: string?
  Algorithm__greedy_group_resolution:
    doc: If this is true, shared peptides will be associated to best proteins only (i.e. become potentially quantifiable razor peptides).
    type: boolean?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: ProteinInference
doc: Protein inference based on an aggregation of the scores of the identified peptides.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - ProteinInference
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json