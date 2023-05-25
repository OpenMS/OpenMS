inputs:
  in:
    doc: "input file "
    type: File
  out:
    doc: output file in mzML format
    type: string
  sampling_rate:
    doc: New sampling rate in m/z dimension (in Th unless ppm flag is set)
    type: double?
  ppm:
    doc: sampling_rate is given in ppm
    type: boolean?
  align_sampling:
    doc: Ensures that sampling is performed equally across the map (same raster is used for all spectra)
    type: boolean?
  min_int_cutoff:
    doc: Intensity cutoff for peaks to be stored in output spectrum (only peaks above this cutoff will be stored, -1 means store all data)
    type: double?
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
label: Resampler
doc: Transforms an LC/MS map into a resampled map or a PNG image.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - Resampler
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json