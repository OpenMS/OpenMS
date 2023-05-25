inputs:
  in_id:
    doc: Peptides with precursor information
    type: File?
  in_text:
    doc: Peptides as text-based file
    type: File?
  in_oligo_params:
    doc: input file with additional model parameters when using the OLIGO kernel
    type: File?
  in_oligo_trainset:
    doc: input file with the used training dataset when using the OLIGO kernel
    type: File?
  svm_model:
    doc: svm model in libsvm format (can be produced by RTModel)
    type: File
  total_gradient_time:
    doc: The time (in seconds) of the gradient (peptide RT prediction)
    type: double?
  max_number_of_peptides:
    doc: The maximum number of peptides considered at once (bigger number will lead to faster results but needs more memory).
    type: long?
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
  out_id__file:
    doc: Output file with peptide RT prediction
    type: string
  out_id__positive:
    doc: Output file in idXML format containing positive predictions for peptide separation prediction - requires 'out_id:negative' to be present as well.
    type: string
  out_id__negative:
    doc: Output file in idXML format containing negative predictions for peptide separation prediction - requires 'out_id:positive' to be present as well.
    type: string
  out_id__rewrite_peptideidentification_rtmz:
    doc: Rewrites each peptideidentification's rt and mz from prediction and calculation (according to the best hit)
    type: boolean?
  out_text__file:
    doc: Output file with predicted RT values
    type: string
outputs:
  out_id__file:
    type: File
    outputBinding:
      glob: $(inputs.out_id__file)
  out_id__positive:
    type: File
    outputBinding:
      glob: $(inputs.out_id__positive)
  out_id__negative:
    type: File
    outputBinding:
      glob: $(inputs.out_id__negative)
  out_text__file:
    type: File
    outputBinding:
      glob: $(inputs.out_text__file)
label: RTPredict
doc: Predicts retention times for peptides using a model trained by RTModel.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - RTPredict
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json