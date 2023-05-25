inputs:
  in:
    doc: "This is the name of the input file (RT prediction). It is assumed that the file type is idXML. Alternatively you can provide a .txt file having a sequence and the corresponding rt per line.\n"
    type: File?
  in_positive:
    doc: "input file with positive examples (peptide separation prediction)\n"
    type: File?
  in_negative:
    doc: "input file with negative examples (peptide separation prediction)\n"
    type: File?
  out:
    doc: "output file: the model in libsvm format"
    type: string
  out_oligo_params:
    doc: output file with additional model parameters when using the OLIGO kernel
    type: string
  out_oligo_trainset:
    doc: output file with the used training dataset when using the OLIGO kernel
    type: string
  svm_type:
    doc: "the type of the svm (NU_SVR or EPSILON_SVR for RT prediction, automatically set\nto C_SVC for separation prediction)\n"
    type: string?
  nu:
    doc: the nu parameter [0..1] of the svm (for nu-SVR)
    type: double?
  p:
    doc: the epsilon parameter of the svm (for epsilon-SVR)
    type: double?
  c:
    doc: the penalty parameter of the svm
    type: double?
  kernel_type:
    doc: the kernel type of the svm
    type: string?
  degree:
    doc: "the degree parameter of the kernel function of the svm (POLY kernel)\n"
    type: long?
  border_length:
    doc: length of the POBK
    type: long?
  max_std:
    doc: max standard deviation for a peptide to be included (if there are several ones for one peptide string)(median is taken)
    type: double?
  k_mer_length:
    doc: k_mer length of the POBK
    type: long?
  sigma:
    doc: sigma of the POBK
    type: double?
  total_gradient_time:
    doc: the time (in seconds) of the gradient (only for RT prediction)
    type: double?
  first_dim_rt:
    doc: if set the model will be built for first_dim_rt
    type: boolean?
  additive_cv:
    doc: "if the step sizes should be interpreted additively (otherwise the actual value is multiplied\nwith the step size to get the new value"
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
    doc: Set to enable Cross-Validation or set to true if the model should just be trained with 1 set of specified parameters.
    type: boolean?
  cv__number_of_runs:
    doc: number of runs for the CV (each run creates a new random partition of the data)
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
  cv__p_start:
    doc: starting point of p
    type: double?
  cv__p_step_size:
    doc: step size point of p
    type: double?
  cv__p_stop:
    doc: stopping point of p
    type: double?
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
label: RTModel
doc: Trains a model for the retention time prediction of peptides from a training set.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - RTModel
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json