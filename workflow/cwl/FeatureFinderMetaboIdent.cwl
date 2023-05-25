inputs:
  in:
    doc: "Input file: LC-MS raw data"
    type: File
  id:
    doc: "Input file: Metabolite identifications"
    type: File
  out:
    doc: "Output file: Features"
    type: string
  lib_out:
    doc: "Output file: Assay library"
    type: string
  chrom_out:
    doc: "Output file: Chromatograms"
    type: string
  trafo_out:
    doc: "Output file: Retention times (expected vs. observed)"
    type: string
  force:
    doc: Overrides tool-specific checks
    type: boolean?
  candidates_out:
    doc: "Optional output file: Feature candidates (before filtering and model fitting)."
    type: string
  debug:
    doc: Sets the debug level
    type: long?
  log:
    doc: Name of log file (created only when specified)
    type: string?
  threads:
    doc: Sets the number of threads allowed to be used by the TOPP tool
    type: long?
  no_progress:
    doc: Disables progress logging to command line
    type: boolean?
  extract__mz_window:
    doc: "m/z window size for chromatogram extraction (unit: ppm if 1 or greater, else Da/Th)"
    type: double?
  extract__rt_window:
    doc: RT window size (in sec.) for chromatogram extraction. If set, this parameter takes precedence over 'extract:rt_quantile'.
    type: double?
  extract__n_isotopes:
    doc: Number of isotopes to include in each peptide assay.
    type: long?
  extract__isotope_pmin:
    doc: Minimum probability for an isotope to be included in the assay for a peptide. If set, this parameter takes precedence over 'extract:n_isotopes'.
    type: double?
  detect__peak_width:
    doc: Expected elution peak width in seconds, for smoothing (Gauss filter). Also determines the RT extration window, unless set explicitly via 'extract:rt_window'.
    type: double?
  detect__min_peak_width:
    doc: Minimum elution peak width. Absolute value in seconds if 1 or greater, else relative to 'peak_width'.
    type: double?
  detect__signal_to_noise:
    doc: Signal-to-noise threshold for OpenSWATH feature detection
    type: double?
  model__type:
    doc: Type of elution model to fit to features
    type: string?
  model__add_zeros:
    doc: Add zero-intensity points outside the feature range to constrain the model fit. This parameter sets the weight given to these points during model fitting; '0' to disable.
    type: double?
  model__unweighted_fit:
    doc: Suppress weighting of mass traces according to theoretical intensities when fitting elution models
    type: boolean?
  model__no_imputation:
    doc: If fitting the elution model fails for a feature, set its intensity to zero instead of imputing a value from the initial intensity estimate
    type: boolean?
  model__each_trace:
    doc: Fit elution model to each individual mass trace
    type: boolean?
  model__check__min_area:
    doc: Lower bound for the area under the curve of a valid elution model
    type: double?
  model__check__boundaries:
    doc: Time points corresponding to this fraction of the elution model height have to be within the data region used for model fitting
    type: double?
  model__check__width:
    doc: Upper limit for acceptable widths of elution models (Gaussian or EGH), expressed in terms of modified (median-based) z-scores. '0' to disable. Not applied to individual mass traces (parameter 'each_trace').
    type: double?
  model__check__asymmetry:
    doc: Upper limit for acceptable asymmetry of elution models (EGH only), expressed in terms of modified (median-based) z-scores. '0' to disable. Not applied to individual mass traces (parameter 'each_trace').
    type: double?
  EMGScoring__max_iteration:
    doc: Maximum number of iterations for EMG fitting.
    type: long?
  EMGScoring__init_mom:
    doc: Alternative initial parameters for fitting through method of moments.
    type: boolean?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  lib_out:
    type: File
    outputBinding:
      glob: $(inputs.lib_out)
  chrom_out:
    type: File
    outputBinding:
      glob: $(inputs.chrom_out)
  trafo_out:
    type: File
    outputBinding:
      glob: $(inputs.trafo_out)
  candidates_out:
    type: File
    outputBinding:
      glob: $(inputs.candidates_out)
label: FeatureFinderMetaboIdent
doc: Detects features in MS1 data based on metabolite identifications.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FeatureFinderMetaboIdent
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json