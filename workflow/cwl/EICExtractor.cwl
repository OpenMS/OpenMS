inputs:
  in:
    doc: Input raw data file
    type: File[]
  in_header:
    doc: "[for Waters data only] Read additional information from _HEADER.TXT. Provide one for each raw input file."
    type: File[]?
  pos:
    doc: Input config file stating where to find signal
    type: File
  rt_tol:
    doc: RT tolerance in [s] for finding max peak (whole RT range around RT middle)
    type: double?
  mz_tol:
    doc: m/z tolerance in [ppm] for finding a peak
    type: double?
  rt_collect:
    doc: "# of scans up & down in RT from highest point for ppm estimation in result"
    type: long?
  out_separator:
    doc: Separator character for output CSV file.
    type: string?
  out:
    doc: Output quantitation file (multiple columns for each input compound)
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
  auto_rt__enabled:
    doc: Automatically detect injection peaks from TIC and quantify all m/z x RT combinations.
    type: boolean?
  auto_rt__FHWM:
    doc: Expected full width at half-maximum of each raw RT peak in [s]. Gaussian smoothing filter with this width is applied to TIC.
    type: double?
  auto_rt__SNThreshold:
    doc: S/N threshold for a smoothed raw peak to pass peak picking. Higher thesholds will result in less peaks.
    type: double?
  auto_rt__out_debug_TIC:
    doc: Optional output file (for first input) containing the smoothed TIC, S/N levels and picked RT positions
    type: string
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  auto_rt__out_debug_TIC:
    type: File
    outputBinding:
      glob: $(inputs.auto_rt__out_debug_TIC)
label: EICExtractor
doc: Extracts intensities from dedicates positions in a LC/MS map
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - EICExtractor
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json