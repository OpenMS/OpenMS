inputs:
  in:
    doc: Mass traces
    type: File
  out:
    doc: output file
    type: string
  min_pearson_correlation:
    doc: Minimal pearson correlation score
    type: double?
  min_peak_nr:
    doc: Minimal peak nr to output pseudo spectra
    type: long?
  max_lag:
    doc: Maximal lag
    type: long?
  max_rt_apex_difference:
    doc: Maximal difference of the apex in retention time
    type: double?
  max_intensity_cutoff:
    doc: Maximal intensity to be added to a spectrum
    type: double?
  add_precursor:
    doc: Add a precursor mass
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
label: ClusterMassTraces
doc: Creates pseudo spectra.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - ClusterMassTraces
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json