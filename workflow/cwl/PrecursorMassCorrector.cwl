inputs:
  in:
    doc: Input mzML file containing the spectra.
    type: File
  out:
    doc: Output mzML file.
    type: string
  feature_in:
    doc: "Input featureXML file, containing features; if set, the MS/MS spectra precursor entries \nwill be matched to the feature m/z values if possible."
    type: File?
  precursor_mass_tolerance:
    doc: "Maximal deviation in Th which is acceptable to be corrected;\nthis value should be set to the instruments selection window."
    type: double?
  max_charge:
    doc: Maximal charge that should be assumed for precursor peaks
    type: long?
  intensity_threshold:
    doc: Intensity threshold value for isotope wavelet feature finder, please look at the documentation of the class for details.
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
label: PrecursorMassCorrector
doc: Corrects the precursor entries of MS/MS spectra, by using MS1 information.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - PrecursorMassCorrector
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json