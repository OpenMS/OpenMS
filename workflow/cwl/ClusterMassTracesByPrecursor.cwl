inputs:
  in_ms1:
    doc: MS1 mass traces
    type: File
  in_swath:
    doc: MS2 / SWATH mass traces
    type: File
  out:
    doc: output file
    type: string
  assign_unassigned_to_all:
    doc: Assign unassigned MS2 fragments to all precursors (only for ms1_centrif)
    type: boolean?
  min_pearson_correlation:
    doc: Minimal pearson correlation score to match elution profiles to each other.
    type: double?
  max_lag:
    doc: Maximal lag (e.g. by how many spectra the peak may be shifted at most). This parameter will depend on your chromatographic setup but a number between 1 and 3 is usually sensible.
    type: long?
  min_nr_ions:
    doc: Minimal number of ions to report a spectrum.
    type: long?
  max_rt_apex_difference:
    doc: Maximal difference of the apex in retention time (in seconds). This is a hard parameter, all profiles further away will not be considered at all.
    type: double?
  swath_lower:
    doc: Swath lower isolation window
    type: double?
  swath_upper:
    doc: Swath upper isolation window
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
label: ClusterMassTracesByPrecursor
doc: Correlate precursor masstraces with fragment ion masstraces in SWATH maps based on their elution profile.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - ClusterMassTracesByPrecursor
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json