inputs:
  in:
    doc: input file
    type: File
  in_type:
    doc: "input file type -- default: determined from file extension or content"
    type: string?
  out:
    doc: Optional output file. If left out, the output is written to the command line.
    type: string
  out_tsv:
    doc: Second optional output file. Tab separated flat text file.
    type: string
  m:
    doc: Show meta information about the whole experiment
    type: boolean?
  p:
    doc: Shows data processing information
    type: boolean?
  s:
    doc: Computes a five-number statistics of intensities, qualities, and widths
    type: boolean?
  d:
    doc: Show detailed listing of all spectra and chromatograms (peak files only)
    type: boolean?
  c:
    doc: Check for corrupt data in the file (peak files only)
    type: boolean?
  v:
    doc: Validate the file only (for mzML, mzData, mzXML, featureXML, idXML, consensusXML, pepXML)
    type: boolean?
  i:
    doc: Check whether a given mzML file contains valid indices (conforming to the indexedmzML standard)
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
  out_tsv:
    type: File
    outputBinding:
      glob: $(inputs.out_tsv)
label: FileInfo
doc: Shows basic information about the file, such as data ranges and file type.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FileInfo
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json