inputs:
  in:
    doc: Input mzML file containing the ER spectra.
    type: File
  pair_in:
    doc: "Pair-file in the format: m/z-light m/z-heavy charge rt"
    type: File
  out:
    doc: Output consensusXML file were the pairs of the feature are written into.
    type: string
  feature_out:
    doc: Output featureXML file, only written if given, skipped otherwise.
    type: string
  precursor_mass_tolerance:
    doc: Precursor mass tolerance which is used for the pair finding and the matching of the given pair m/z values to the features.
    type: double?
  RT_tolerance:
    doc: Maximal deviation in RT dimension in seconds a feature can have when comparing to the RT values given in the pair file
    type: double?
  max_charge:
    doc: Maximal charge state features should be search for.
    type: long?
  intensity_threshold:
    doc: Intensity threshold, for the meaning see the documentation of the IsotopeWaveletFeatureFinder documentation.
    type: double?
  max_isotope:
    doc: Max isotope of the isotope distribution to be considered
    type: long?
  expansion_range:
    doc: The range that is used to extend the isotope distribution with null intensity peaks in Th.
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
  feature_out:
    type: File
    outputBinding:
      glob: $(inputs.feature_out)
label: ERPairFinder
doc: Util which can be used to evaluate pair ratios on enhanced resolution (zoom) scans.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - ERPairFinder
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json