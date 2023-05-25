inputs:
  in:
    doc: "input raw data file "
    type: File
  out:
    doc: "output raw data file "
    type: string
  processOption:
    doc: Whether to load all data and process them in-memory or whether to process the data on the fly (lowmemory) without loading the whole file into memory first
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
  algorithm__gaussian_width:
    doc: Use a gaussian filter width which has approximately the same width as your mass peaks (FWHM in m/z).
    type: double?
  algorithm__ppm_tolerance:
    doc: "Gaussian width, depending on the m/z position.\nThe higher the value, the wider the peak and therefore the wider the gaussian."
    type: double?
  algorithm__use_ppm_tolerance:
    doc: If true, instead of the gaussian_width value, the ppm_tolerance is used. The gaussian is calculated in each step anew, so this is much slower.
    type: boolean?
  algorithm__write_log_messages:
    doc: "true: Warn if no signal was found by the Gauss filter algorithm."
    type: boolean?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: NoiseFilterGaussian
doc: Removes noise from profile spectra by using Gaussian filter (on uniform as well as non-uniform data).
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - NoiseFilterGaussian
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json