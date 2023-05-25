inputs:
  in:
    doc: Input files separated by blank
    type: File[]
  tr:
    doc: transition file ('TraML' or 'csv')
    type: File
  rt_norm:
    doc: RT normalization file (how to map the RTs of this run to the ones stored in the library)
    type: File?
  out:
    doc: output file
    type: string
  min_upper_edge_dist:
    doc: Minimal distance to the edge to still consider a precursor, in Thomson
    type: double?
  rt_window:
    doc: Extraction window in RT dimension (-1 means extract over the whole range). This is the full window size, e.g. a value of 1000 seconds would extract 500 seconds on either side.
    type: double?
  ion_mobility_window:
    doc: Extraction window in ion mobility dimension (in milliseconds). This is the full window size, e.g. a value of 10 milliseconds would extract 5 milliseconds on either side.
    type: double?
  mz_window:
    doc: Extraction window in m/z dimension (in Thomson, to use ppm see -ppm flag). This is the full window size, e.g. 100 ppm would extract 50 ppm on either side.
    type: double?
  ppm:
    doc: m/z extraction_window is in ppm
    type: boolean?
  is_swath:
    doc: Set this flag if the data is SWATH data
    type: boolean?
  extract_MS1:
    doc: Extract the MS1 transitions based on the precursor values in the TraML file (useful for extracting MS1 XIC)
    type: boolean?
  extraction_function:
    doc: Function used to extract the signal
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
  model__type:
    doc: Type of model
    type: string?
  model__symmetric_regression:
    doc: "Only for 'linear' model: Perform linear regression on 'y - x' vs. 'y + x', instead of on 'y' vs. 'x'."
    type: boolean?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: OpenSwathChromatogramExtractor
doc: Extract chromatograms (XIC) from a MS2 map file.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - OpenSwathChromatogramExtractor
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json