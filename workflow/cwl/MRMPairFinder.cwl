inputs:
  in:
    doc: Input featureXML file containing the features of the MRM experiment spectra.
    type: File
  pair_in:
    doc: "Pair-file in the format: prec-m/z-light prec-m/z-heavy frag-m/z-light frag-m/z-heavy rt"
    type: File
  out:
    doc: Output consensusXML file were the pairs of the features will be written to.
    type: string
  feature_out:
    doc: Output featureXML file, only written if given, skipped otherwise.
    type: string
  mass_tolerance:
    doc: Precursor mass tolerance which is used for the pair finding and the matching of the given pair m/z values to the features.
    type: double?
  RT_tolerance:
    doc: Maximal deviation in RT dimension in seconds a feature can have when comparing to the RT values given in the pair file.
    type: double?
  RT_pair_tolerance:
    doc: Maximal deviation in RT dimension in seconds the two partners of a pair is allowed to have.
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
label: MRMPairFinder
doc: Util which can be used to evaluate labeled pair ratios on MRM features.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MRMPairFinder
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json