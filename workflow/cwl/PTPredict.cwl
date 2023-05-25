inputs:
  in:
    doc: "input file "
    type: File
  in_oligo_params:
    doc: input file with additional model parameters when using the OLIGO kernel
    type: File?
  in_oligo_trainset:
    doc: input file with the used training dataset when using the OLIGO kernel
    type: File?
  out:
    doc: "output file\n"
    type: string
  svm_model:
    doc: svm model in libsvm format (can be produced by PTModel)
    type: File
  max_number_of_peptides:
    doc: "the maximum number of peptides considered at once (bigger number will lead to faster results but needs more memory).\n"
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
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: PTPredict
doc: predicts the likelihood of peptides to be proteotypic via svm_model which is trained by PTModel
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - PTPredict
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json