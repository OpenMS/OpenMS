inputs:
  in:
    doc: input file
    type: File
  out:
    doc: output file
    type: string
  algorithm_type:
    doc: "The normalization algorithm that is applied. 'robust_regression' scales each map by a fator computed from the ratios of non-differential background features (as determined by the ratio_threshold parameter), 'quantile' performs quantile normalization, 'median' scales all maps to the same median intensity, 'median_shift' shifts the median instead of scaling (WARNING: if you have regular, log-normal MS data, 'median_shift' is probably the wrong choice. Use only if you know what you're doing!)"
    type: string?
  ratio_threshold:
    doc: "Only for 'robust_regression': the parameter is used to distinguish between non-outliers (ratio_threshold < intensity ratio < 1/ratio_threshold) and outliers."
    type: double?
  accession_filter:
    doc: Use only features with accessions (partially) matching this regular expression for computing the normalization factors. Useful, e.g., if you have known house keeping proteins in your samples. When this parameter is empty or the regular expression matches the empty string, all features are used (even those without an ID). No effect if quantile normalization is used.
    type: string?
  description_filter:
    doc: Use only features with description (partially) matching this regular expression for computing the normalization factors. Useful, e.g., if you have known house keeping proteins in your samples. When this parameter is empty or the regular expression matches the empty string, all features are used (even those without an ID). No effect if quantile normalization is used.
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
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: ConsensusMapNormalizer
doc: Normalizes maps of one consensusXML file
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - ConsensusMapNormalizer
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json