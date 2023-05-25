inputs:
  in_positive:
    doc: input file with positive examples
    type: File
  in_negative:
    doc: input file with negative examples
    type: File
  out:
    doc: "output file: the model in libsvm format"
    type: string
  out_oligo_params:
    doc: output file with additional model parameters when using the OLIGO kernel
    type: string
  out_oligo_trainset:
    doc: output file with the used training dataset when using the OLIGO kernel
    type: string
  c:
    doc: the penalty parameter of the svm
    type: double?
  svm_type:
    doc: the type of the svm (NU_SVC or C_SVC)
    type: string?
  nu:
    doc: the nu parameter [0..1] of the svm (for nu-SVR)
    type: double?
  kernel_type:
    doc: the kernel type of the svm
    type: string?
  degree:
    doc: the degree parameter of the kernel function of the svm (POLY kernel)
    type: long?
  border_length:
    doc: length of the POBK
    type: long?
  k_mer_length:
    doc: k_mer length of the POBK
    type: long?
  sigma:
    doc: sigma of the POBK
    type: double?
  max_positive_count:
    doc: quantity of positive samples for training (randomly chosen if smaller than available quantity)
    type: long?
  max_negative_count:
    doc: quantity of positive samples for training (randomly chosen if smaller than available quantity)
    type: long?
  redundant:
    doc: if the input sets are redundant and the redundant peptides should occur more than once in the training set, this flag has to be set
    type: boolean?
  additive_cv:
    doc: if the step sizes should be interpreted additively (otherwise the actual value is multiplied with the step size to get the new value
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
  cv__skip_cv:
    doc: Has to be set if the cv should be skipped and the model should just be trained with the specified parameters.
    type: boolean?
  cv__number_of_runs:
    doc: number of runs for the CV
    type: long?
  cv__number_of_partitions:
    doc: number of CV partitions
    type: long?
  cv__degree_start:
    doc: starting point of degree
    type: long?
  cv__degree_step_size:
    doc: step size point of degree
    type: long?
  cv__degree_stop:
    doc: stopping point of degree
    type: long?
  cv__c_start:
    doc: starting point of c
    type: double?
  cv__c_step_size:
    doc: step size of c
    type: double?
  cv__c_stop:
    doc: stopping point of c
    type: double?
  cv__nu_start:
    doc: starting point of nu
    type: double?
  cv__nu_step_size:
    doc: step size of nu
    type: double?
  cv__nu_stop:
    doc: stopping point of nu
    type: double?
  cv__sigma_start:
    doc: starting point of sigma
    type: double?
  cv__sigma_step_size:
    doc: step size of sigma
    type: double?
  cv__sigma_stop:
    doc: stopping point of sigma
    type: double?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_oligo_params:
    type: File
    outputBinding:
      glob: $(inputs.out_oligo_params)
  out_oligo_trainset:
    type: File
    outputBinding:
      glob: $(inputs.out_oligo_trainset)
label: PTModel
doc: Trains a model for the prediction of proteotypic peptides from a training set.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - PTModel
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json