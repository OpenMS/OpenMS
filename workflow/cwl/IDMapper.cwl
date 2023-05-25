inputs:
  id:
    doc: Protein/peptide identifications file
    type: File
  in:
    doc: Feature map/consensus map file
    type: File
  out:
    doc: Output file (the format depends on the input file format).
    type: string
  rt_tolerance:
    doc: "RT tolerance (in seconds) for the matching of peptide identifications and (consensus) features.\nTolerance is understood as 'plus or minus x', so the matching range increases by twice the given value."
    type: double?
  mz_tolerance:
    doc: "m/z tolerance (in ppm or Da) for the matching of peptide identifications and (consensus) features.\nTolerance is understood as 'plus or minus x', so the matching range increases by twice the given value."
    type: double?
  mz_measure:
    doc: Unit of 'mz_tolerance'.
    type: string?
  mz_reference:
    doc: "Source of m/z values for peptide identifications. If 'precursor', the precursor-m/z from the idXML is used. If 'peptide',\nmasses are computed from the sequences of peptide hits; in this case, an identification matches if any of its hits matches.\n('peptide' should be used together with 'feature:use_centroid_mz' to avoid false-positive matches.)"
    type: string?
  ignore_charge:
    doc: "For feature/consensus maps: Assign an ID independently of whether its charge state matches that of the (consensus) feature."
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
  feature__use_centroid_rt:
    doc: Use the RT coordinates of the feature centroids for matching, instead of the RT ranges of the features/mass traces.
    type: boolean?
  feature__use_centroid_mz:
    doc: "Use the m/z coordinates of the feature centroids for matching, instead of the m/z ranges of the features/mass traces.\n(If you choose 'peptide' as 'mz_reference', you should usually set this flag to avoid false-positive matches.)"
    type: string?
  consensus__use_subelements:
    doc: Match using RT and m/z of sub-features instead of consensus RT and m/z. A consensus feature matches if any of its sub-features matches.
    type: boolean?
  consensus__annotate_ids_with_subelements:
    doc: Store the map index of the sub-feature in the peptide ID.
    type: boolean?
  spectra__in:
    doc: MS run used to annotated unidentified spectra to features or consensus features.
    type: File?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: IDMapper
doc: Assigns protein/peptide identifications to features or consensus features.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - IDMapper
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json