inputs:
  control:
    doc: input mzML file
    type: File
  treatment:
    doc: input mzML file
    type: File
  fold_change:
    doc: fold change between XICs
    type: double?
  rt_tol:
    doc: RT tolerance in [s] for finding max peak (whole RT range around RT middle)
    type: double?
  mz_tol:
    doc: m/z tolerance in [ppm] for finding a peak
    type: double?
  out:
    doc: output of the treatment file after XIC filtering.
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
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: RNPxlXICFilter
doc: Remove MS2 spectra from treatment based on the fold change between control and treatment.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - RNPxlXICFilter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json