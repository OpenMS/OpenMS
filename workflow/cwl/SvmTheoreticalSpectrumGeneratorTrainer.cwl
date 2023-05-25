inputs:
  in_spectra:
    doc: Input Training Spectra in mzML
    type: File
  in_identifications:
    doc: Input file with corresponding sequences in idXML
    type: File
  model_output_file:
    doc: Name for output files. For each ion_type one file <filename>_residue_loss_charge.svm and one <filename>.info which has to be passed to the SvmTheoretical SpectrumGenerator
    type: string
  precursor_charge:
    doc: Precursor charge state used for model training
    type: long?
  write_training_files:
    doc: No models are trained but input training files for libSVM command line tools are produced
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
  algorithm__number_intensity_levels:
    doc: The number of intensity bins (for secondary type models)
    type: long?
  algorithm__number_regions:
    doc: The number of regions each spectrum is split to (for secondary type models)
    type: long?
  algorithm__parent_tolerance:
    doc: The maximum difference between theoretical and experimental parent mass to accept training spectrum
    type: double?
  algorithm__peak_tolerance:
    doc: The maximum mass error for a peak to the expected mass of some ion type
    type: double?
  algorithm__add_b_ions:
    doc: Train simulator for b-ions
    type: string?
  algorithm__add_y_ions:
    doc: Train simulator for y-ions
    type: string?
  algorithm__add_a_ions:
    doc: Train simulator for a-ions
    type: boolean?
  algorithm__add_c_ions:
    doc: Train simulator for c-ions
    type: boolean?
  algorithm__add_x_ions:
    doc: Train simulator for x-ions
    type: boolean?
  algorithm__add_z_ions:
    doc: Train simulator for z-ions
    type: boolean?
  algorithm__add_losses:
    doc: Train simulator for neutral losses of H2O and NH3 for b-ions and y-ions
    type: boolean?
  algorithm__add_b2_ions:
    doc: Train simulator for doubly charged b-ions
    type: boolean?
  algorithm__add_y2_ions:
    doc: Train simulator for double charged y-ions
    type: boolean?
  algorithm__svm__svc_type:
    doc: "Type of the SVC: 0=C_SVC 1=NU_SVC"
    type: long?
  algorithm__svm__svr_type:
    doc: "Type of the SVR: 0=EPSILON_SVR 1=NU_SVR"
    type: long?
  algorithm__svm__scaling:
    doc: Apply scaling of feature values
    type: string?
  algorithm__svm__scaling_lower:
    doc: Lower bound for scaling
    type: double?
  algorithm__svm__scaling_upper:
    doc: Upper bound for scaling
    type: double?
  algorithm__svm__n_fold:
    doc: n_fold cross validation is performed
    type: long?
  algorithm__svm__grid:
    doc: Perform grid search
    type: boolean?
  algorithm__svm__additive_cv:
    doc: Additive step size (if false multiplicative)
    type: boolean?
  algorithm__svm__svc__kernel_type:
    doc: "Type of the kernel:  0=LINEAR 1=POLY 2=RBF 3=SIGMOID"
    type: long?
  algorithm__svm__svc__degree:
    doc: For POLY
    type: long?
  algorithm__svm__svc__gamma:
    doc: For POLY/RBF/SIGMOID
    type: double?
  algorithm__svm__svc__C:
    doc: Cost of constraint violation
    type: double?
  algorithm__svm__svc__nu:
    doc: For NU_SVC, ONE_CLASS and NU_SVR
    type: double?
  algorithm__svm__svc__balancing:
    doc: Use class balanced SVC training
    type: string?
  algorithm__svm__svc__degree_start:
    doc: starting point of degree
    type: long?
  algorithm__svm__svc__degree_step_size:
    doc: step size point of degree
    type: long?
  algorithm__svm__svc__degree_stop:
    doc: stopping point of degree
    type: long?
  algorithm__svm__svc__gamma_start:
    doc: starting point of gamma
    type: double?
  algorithm__svm__svc__gamma_step_size:
    doc: step size point of gamma
    type: long?
  algorithm__svm__svc__gamma_stop:
    doc: stopping point of gamma
    type: double?
  algorithm__svm__svc__c_start:
    doc: starting point of c
    type: double?
  algorithm__svm__svc__c_step_size:
    doc: step size of c
    type: long?
  algorithm__svm__svc__c_stop:
    doc: stopping point of c
    type: long?
  algorithm__svm__svc__nu_start:
    doc: starting point of nu
    type: double?
  algorithm__svm__svc__nu_step_size:
    doc: step size of nu
    type: long?
  algorithm__svm__svc__nu_stop:
    doc: stopping point of nu
    type: double?
  algorithm__svm__svr__kernel_type:
    doc: "Type of the kernel:  0=LINEAR 1=POLY 2=RBF 3=SIGMOID"
    type: long?
  algorithm__svm__svr__degree:
    doc: For POLY
    type: long?
  algorithm__svm__svr__gamma:
    doc: For POLY/RBF/SIGMOID
    type: double?
  algorithm__svm__svr__C:
    doc: Cost of constraint violation
    type: double?
  algorithm__svm__svr__p:
    doc: The epsilon for the loss function in epsilon-SVR
    type: double?
  algorithm__svm__svr__nu:
    doc: For NU_SVC, ONE_CLASS and NU_SVR
    type: double?
  algorithm__svm__svr__degree_start:
    doc: starting point of degree
    type: long?
  algorithm__svm__svr__degree_step_size:
    doc: step size point of degree
    type: long?
  algorithm__svm__svr__degree_stop:
    doc: stopping point of degree
    type: long?
  algorithm__svm__svr__gamma_start:
    doc: starting point of gamma
    type: double?
  algorithm__svm__svr__gamma_step_size:
    doc: step size point of gamma
    type: long?
  algorithm__svm__svr__gamma_stop:
    doc: stopping point of gamma
    type: double?
  algorithm__svm__svr__p_start:
    doc: starting point of p
    type: double?
  algorithm__svm__svr__p_step_size:
    doc: step size point of p
    type: long?
  algorithm__svm__svr__p_stop:
    doc: stopping point of p
    type: double?
  algorithm__svm__svr__c_start:
    doc: starting point of c
    type: double?
  algorithm__svm__svr__c_step_size:
    doc: step size of c
    type: long?
  algorithm__svm__svr__c_stop:
    doc: stopping point of c
    type: long?
  algorithm__svm__svr__nu_start:
    doc: starting point of nu
    type: double?
  algorithm__svm__svr__nu_step_size:
    doc: step size of nu
    type: long?
  algorithm__svm__svr__nu_stop:
    doc: stopping point of nu
    type: double?
outputs:
  model_output_file:
    type: File
    outputBinding:
      glob: $(inputs.model_output_file)
label: SvmTheoreticalSpectrumGeneratorTrainer
doc: Trainer for SVM models as input for SvmTheoreticalSpectrumGenerator
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - SvmTheoreticalSpectrumGeneratorTrainer
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json