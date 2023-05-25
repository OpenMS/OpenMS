inputs:
  in:
    doc: input file
    type: File
  out:
    doc: output file
    type: string
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
  algorithm__min_rt_distance:
    doc: Minimal distance of MRM features in seconds.
    type: double?
  algorithm__min_num_peaks_per_feature:
    doc: Minimal number of peaks which are needed for a single feature
    type: long?
  algorithm__min_signal_to_noise_ratio:
    doc: Minimal S/N ratio a peak must have to be taken into account. Set to zero if the MRM-traces contains mostly signals, and no noise.
    type: double?
  algorithm__write_debug_files:
    doc: If set to true, for each feature a plot will be created, in the subdirectory 'debug'
    type: boolean?
  algorithm__resample_traces:
    doc: If set to true, each trace, which is in this case a part of the MRM monitoring trace with signal is resampled, using the minimal distance of two data points in RT dimension
    type: boolean?
  algorithm__write_debuginfo:
    doc: If set to true, debug messages are written, the output can be somewhat lengthy.
    type: boolean?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: FeatureFinderMRM
doc: Detects two-dimensional features in LC-MS data.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FeatureFinderMRM
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json