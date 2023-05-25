inputs:
  in:
    doc: Input file (mzML)
    type: File
  out:
    doc: Default output tsv file containing deconvolved features
    type: string
  out_mzml:
    doc: Output mzml file containing deconvolved spectra (of all MS levels)
    type: string
  out_annotated_mzml:
    doc: Output mzml file containing annotated spectra. For each annotated peak, monoisotopic mass, charge, and isotope index are stored as meta data. Unannotated peaks are also copied as well without meta data.
    type: string
  out_promex:
    doc: Output ms1ft (promex compatible) file containing deconvolved spectra. Only for MS1 level
    type: string
  min_precursor_snr:
    doc: Minimum precursor SNR (SNR within the precursor envelope range) for identification. Similar to precursor interference level, but more stringent.When FLASHIda log file is used, this parameter is ignored. Applied only for topFD msalign outputs.
    type: double?
  mzml_mass_charge:
    doc: Charge status of deconvolved masses in mzml output (specified by out_mzml)
    type: long?
  preceding_MS1_count:
    doc: Specifies the number of preceding MS1 spectra for MS2 precursor determination. In TDP, the precursor peak of a MS2 spectrum may not belong to any deconvolved masses in the MS1 spectrum immediately preceding the MS2 spectrum. Increasing this parameter to N allows for the search for the deconvolved masses in the N preceding MS1 spectra from the MS2 spectrum, increasing the chance that its precursor is deconvolved.
    type: long?
  write_detail:
    doc: To write peak information per deconvolved mass in detail or not in tsv files for deconvolved spectra. If set to 1, all peak information (m/z, intensity, charge and isotope index) per mass is reported.
    type: long?
  max_MS_level:
    doc: Maximum MS level (inclusive) for deconvolution.
    type: long?
  forced_MS_level:
    doc: If set to an integer N, MS level of all spectra will be set to N regardless of original MS level. Useful when deconvolving datasets containing only MS2 spectra.
    type: long?
  merging_method:
    doc: "Method for spectra merging before deconvolution. 0: No merging 1: Average gaussian method to perform moving gaussian averaging of spectra per MS level. Effective to increase proteoform ID sensitivity (in particular for Q-TOF datasets). 2: Block method to perform merging of all spectra into a single one per MS level (e.g., for NativeMS datasets)"
    type: long?
  report_FDR:
    doc: Report qvalues (roughly, mass-wise FDR) for deconvolved masses in the tsv files for deconvolved spectra. Decoy masses to calculate qvalues and FDR are also reported. Beta version.
    type: long?
  use_RNA_averagine:
    doc: If set to 1, RNA averagine model is used
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
  Algorithm__tol:
    doc: "ppm tolerance for MS1, MS2, ... "
    type: double[]?
  Algorithm__min_mass:
    doc: Minimum mass (Da)
    type: double?
  Algorithm__max_mass:
    doc: Maximum mass (Da)
    type: double?
  Algorithm__min_charge:
    doc: Minimum charge state for MS1 spectra (can be negative for negative mode)
    type: long?
  Algorithm__max_charge:
    doc: Maximum charge state for MS1 spectra (can be negative for negative mode)
    type: long?
  Algorithm__min_mz:
    doc: If set to positive value, minimum m/z to deconvolve.
    type: double?
  Algorithm__max_mz:
    doc: If set to positive value, maximum m/z to deconvolve.
    type: double?
  Algorithm__min_rt:
    doc: If set to positive value, minimum RT to deconvolve.
    type: double?
  Algorithm__max_rt:
    doc: If set to positive value, maximum RT to deconvolve.
    type: double?
  Algorithm__isolation_window:
    doc: Default isolation window with. If the input mzML file does not contain isolation window width information, this width will be used.
    type: double?
  Algorithm__min_isotope_cosine:
    doc: Cosine similarity thresholds between avg. and observed isotope patterns for MS1, 2, ... (e.g., -min_isotope_cosine 0.8 0.6 to specify 0.8 and 0.6 for MS1 and MS2, respectively)
    type: double[]?
  Algorithm__allowed_isotope_error:
    doc: Allowed isotope index error for decoy and qvalue report. If it is set to 1, for example, +-1 isotope errors are not counted as false. Beta version.
    type: long?
  Algorithm__min_intensity:
    doc: Intensity threshold
    type: double?
  FeatureTracing__mass_error_ppm:
    doc: Feature tracing mass ppm tolerance. When negative, MS1 tolerance for mass deconvolution will be used (e.g., 16 ppm is used when -Algorithm:tol 16).
    type: double?
  FeatureTracing__quant_method:
    doc: Method of quantification for mass traces. For LC data 'area' is recommended, 'median' for direct injection data. 'max_height' simply uses the most intense peak in the trace.
    type: string?
  FeatureTracing__min_sample_rate:
    doc: Minimum fraction of scans along the feature trace that must contain a peak. To raise feature detection sensitivity, lower this value close to 0.
    type: double?
  FeatureTracing__min_trace_length:
    doc: Minimum expected length of a mass trace (in seconds).
    type: double?
  FeatureTracing__max_trace_length:
    doc: Maximum expected length of a mass trace (in seconds). Set to a negative value to disable maximal length check during mass trace detection.
    type: double?
  FeatureTracing__min_isotope_cosine:
    doc: Cosine similarity threshold between avg. and observed isotope pattern for mass features. if not set, controlled by -Algorithm:min_isotope_cosine_ option
    type: double?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_mzml:
    type: File
    outputBinding:
      glob: $(inputs.out_mzml)
  out_annotated_mzml:
    type: File
    outputBinding:
      glob: $(inputs.out_annotated_mzml)
  out_promex:
    type: File
    outputBinding:
      glob: $(inputs.out_promex)
label: FLASHDeconv
doc: Ultra-fast high-quality deconvolution enables online processing of top-down MS data
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FLASHDeconv
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json