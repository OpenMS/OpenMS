inputs:
  in:
    doc: Input mzML file
    type: File
  in_type:
    doc: "input file type -- default: determined from file extension or content\n"
    type: string?
  out:
    doc: Output file
    type: string
  out_type:
    doc: "Output file type -- default: determined from file extension or content\nNote: that not all conversion paths work or make sense."
    type: string?
  convert_back:
    doc: Convert back to mzML
    type: boolean?
  lossy_compression:
    doc: "Use numpress compression to achieve optimally small file size (attention: may cause small loss of precision; only for mzML data)."
    type: string?
  full_meta:
    doc: Write full meta information into sqMass file (may require large amounts of memory)
    type: string?
  lossy_mass_accuracy:
    doc: Desired (absolute) m/z accuracy for lossy compression (e.g. use 0.0001 for a mass accuracy of 0.2 ppm at 500 m/z, default uses -1.0 for maximal accuracy).
    type: double?
  process_lowmemory:
    doc: "Whether to process the file on the fly without loading the whole file into memory first (only for conversions of mzXML/mzML to mzML).\nNote: this flag will prevent conversion from spectra to chromatograms."
    type: boolean?
  lowmem_batchsize:
    doc: The batch size of the low memory conversion
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
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: OpenSwathMzMLFileCacher
doc: This tool caches the spectra and chromatogram data of an mzML to disk.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - OpenSwathMzMLFileCacher
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json