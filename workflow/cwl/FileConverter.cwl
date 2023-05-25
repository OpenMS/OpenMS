inputs:
  in:
    doc: Input file to convert.
    type: File
  in_type:
    doc: "Input file type -- default: determined from file extension or content\n"
    type: string?
  UID_postprocessing:
    doc: "unique ID post-processing for output data.\n'none' keeps current IDs even if invalid.\n'ensure' keeps current IDs but reassigns invalid ones.\n'reassign' assigns new unique IDs."
    type: string?
  out:
    doc: Output file
    type: string
  out_type:
    doc: "Output file type -- default: determined from file extension or content\nNote: that not all conversion paths work or make sense."
    type: string?
  TIC_DTA2D:
    doc: Export the TIC instead of the entire experiment in mzML/mzData/mzXML -> DTA2D conversions.
    type: boolean?
  MGF_compact:
    doc: Use a more compact format when writing MGF (no zero-intensity peaks, limited number of decimal places)
    type: boolean?
  force_MaxQuant_compatibility:
    doc: "[mzXML output only] Make sure that MaxQuant can read the mzXML and set the msManufacturer to 'Thermo Scientific'."
    type: boolean?
  force_TPP_compatibility:
    doc: "[mzML output only] Make sure that TPP parsers can read the mzML and the precursor ion m/z in the file (otherwise it will be set to zero by the TPP)."
    type: boolean?
  convert_to_chromatograms:
    doc: "[mzML output only] Assumes that the provided spectra represent data in SRM mode or targeted MS1 mode and converts them to chromatogram data."
    type: boolean?
  change_im_format:
    doc: "[mzML output only] How to store ion mobility scans (none: no change in format; multiple_spectra: store each IM frame as multiple scans (one per drift time value); concatenated: store whole frame as single scan with IM values in a FloatDataArray"
    type: string?
  write_scan_index:
    doc: Append an index when writing mzML or mzXML files. Some external tools might rely on it.
    type: string?
  lossy_compression:
    doc: "Use numpress compression to achieve optimally small file size using linear compression for m/z domain and slof for intensity and float data arrays (attention: may cause small loss of precision; only for mzML data)."
    type: boolean?
  lossy_mass_accuracy:
    doc: Desired (absolute) m/z accuracy for lossy compression (e.g. use 0.0001 for a mass accuracy of 0.2 ppm at 500 m/z, default uses -1.0 for maximal accuracy).
    type: double?
  process_lowmemory:
    doc: "Whether to process the file on the fly without loading the whole file into memory first (only for conversions of mzXML/mzML to mzML).\nNote: this flag will prevent conversion from spectra to chromatograms."
    type: boolean?
  NET_executable:
    doc: The .NET framework executable. Only required on linux and mac.
    type: File?
  ThermoRaw_executable:
    doc: The ThermoRawFileParser executable.
    type: File?
  no_peak_picking:
    doc: Disables vendor peak picking for raw files.
    type: boolean?
  no_zlib_compression:
    doc: Disables zlib compression for raw file conversion. Enables compatibility with some tools that do not support compressed input files, e.g. X!Tandem.
    type: boolean?
  include_noise:
    doc: Include noise data in mzML output.
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
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: FileConverter
doc: Converts between different MS file formats.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FileConverter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json