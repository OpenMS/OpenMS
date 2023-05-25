inputs:
  in:
    doc: Input mzML file list, containing the spectra.
    type: File[]
  id_in:
    doc: Input idXML file list, containing the identifications.
    type: File[]
  out_precursor:
    doc: Output file which contains the deviations from the precursors
    type: string
  precursor_error_ppm:
    doc: If this flag is used, the precursor mass tolerances are estimated in ppm instead of Da.
    type: boolean?
  out_fragment:
    doc: Output file which contains the fragment ion m/z deviations
    type: string
  fragment_error_ppm:
    doc: If this flag is used, the fragment mass tolerances are estimated in ppm instead of Da.
    type: boolean?
  fragment_mass_tolerance:
    doc: Maximal fragment mass tolerance which is allowed for MS/MS spectra, used for the calculation of matching ions.
    type: double?
  number_of_bins:
    doc: Number of bins that should be used to calculate the histograms for the fitting.
    type: long?
  out_precursor_fit:
    doc: Gaussian fit to the histogram of mass deviations from the precursors.
    type: string
  out_fragment_fit:
    doc: Gaussian fit to the histogram of mass deviations from the fragments.
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
  out_precursor:
    type: File
    outputBinding:
      glob: $(inputs.out_precursor)
  out_fragment:
    type: File
    outputBinding:
      glob: $(inputs.out_fragment)
  out_precursor_fit:
    type: File
    outputBinding:
      glob: $(inputs.out_precursor_fit)
  out_fragment_fit:
    type: File
    outputBinding:
      glob: $(inputs.out_fragment_fit)
label: IDMassAccuracy
doc: Calculates a distribution of the mass error from given mass spectra and IDs.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - IDMassAccuracy
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json