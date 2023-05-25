inputs:
  in:
    doc: "input file "
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
  algorithm__signal_to_noise_:
    doc: Signal to noise value, each peak is required to be above this value (turn off by setting it to 0.0)
    type: double?
  algorithm__peak_width:
    doc: Expected peak width half width in Dalton - peaks will be extended until this half width is reached (even if the intensitity is increasing). In conjunction with check_width_internally it will also be used to remove peaks whose spacing is larger than this value.
    type: double?
  algorithm__spacing_difference:
    doc: Difference between peaks in multiples of the minimal difference to continue. The higher this value is set, the further apart peaks are allowed to be to still extend a peak. E.g. if the value is set to 1.5 and in a current peak the minimal spacing between peaks is 10 mDa, then only peaks at most 15 mDa apart will be added to the peak.
    type: double?
  algorithm__sn_bin_count_:
    doc: Bin count for the Signal to Noise estimation.
    type: long?
  algorithm__nr_iterations_:
    doc: Nr of iterations to perform (how many times the peaks are re-centered).
    type: long?
  algorithm__sn_win_len_:
    doc: Window length for the Signal to Noise estimation.
    type: double?
  algorithm__check_width_internally:
    doc: Delete peaks where the spacing is larger than the peak width (should be set to true to avoid artefacts)
    type: boolean?
  algorithm__ms1_only:
    doc: Only do MS1
    type: boolean?
  algorithm__clear_meta_data:
    doc: Delete meta data about peak width
    type: boolean?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: PeakPickerIterative
doc: Finds mass spectrometric peaks in profile mass spectra.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - PeakPickerIterative
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json