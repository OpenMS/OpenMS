inputs:
  in:
    doc: Input file
    type: File
  out:
    doc: "Prefix for output files ('_part1of2.mzML' etc. will be appended; default: same as 'in' without the file extension)"
    type: string
  parts:
    doc: Number of parts to split into (takes precedence over 'size' if set)
    type: long?
  size:
    doc: Approximate upper limit for resulting file sizes (in 'unit')
    type: long?
  unit:
    doc: Unit for 'size' (base 1024)
    type: string?
  no_chrom:
    doc: Remove chromatograms, keep only spectra.
    type: boolean?
  no_spec:
    doc: Remove spectra, keep only chromatograms.
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
      glob: $(inputs.out)*
label: MzMLSplitter
doc: Splits an mzML file into multiple parts
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MzMLSplitter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json