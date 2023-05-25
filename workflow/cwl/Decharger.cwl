inputs:
  in:
    doc: "input file "
    type: File
  out_cm:
    doc: output consensus map
    type: string
  out_fm:
    doc: output feature map
    type: string
  outpairs:
    doc: output file
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
  algorithm__FeatureDeconvolution__charge_min:
    doc: Minimal possible charge
    type: long?
  algorithm__FeatureDeconvolution__charge_max:
    doc: Maximal possible charge
    type: long?
  algorithm__FeatureDeconvolution__charge_span_max:
    doc: Maximal range of charges for a single analyte, i.e. observing q1=[5,6,7] implies span=3. Setting this to 1 will only find adduct variants of the same charge
    type: long?
  algorithm__FeatureDeconvolution__q_try:
    doc: Try different values of charge for each feature according to the above settings ('heuristic' [does not test all charges, just the likely ones] or 'all' ), or leave feature charge untouched ('feature').
    type: string?
  algorithm__FeatureDeconvolution__retention_max_diff:
    doc: Maximum allowed RT difference between any two features if their relation shall be determined
    type: double?
  algorithm__FeatureDeconvolution__retention_max_diff_local:
    doc: Maximum allowed RT difference between between two co-features, after adduct shifts have been accounted for (if you do not have any adduct shifts, this value should be equal to 'retention_max_diff', otherwise it should be smaller!)
    type: double?
  algorithm__FeatureDeconvolution__mass_max_diff:
    doc: Maximum allowed mass difference [in Th] for a single feature.
    type: double?
  algorithm__FeatureDeconvolution__potential_adducts:
    doc: "Adducts used to explain mass differences in format: 'Element:Charge(+/-):Probability[:RTShift[:Label]]', i.e. the number of '+' or '-' indicate the charge, e.g. 'Ca:++:0.5' indicates +2. Probabilites have to be in (0,1]. RTShift param is optional and indicates the expected RT shift caused by this adduct, e.g. '(2)H4H-4:0:1:-3' indicates a 4 deuterium label, which causes early elution by 3 seconds. As a fifth parameter you can add a label which is tagged on every feature which has this adduct. This also determines the map number in the consensus file."
    type: string[]?
  algorithm__FeatureDeconvolution__max_neutrals:
    doc: Maximal number of neutral adducts(q=0) allowed. Add them in the 'potential_adducts' section!
    type: long?
  algorithm__FeatureDeconvolution__max_minority_bound:
    doc: Maximum count of the least probable adduct (according to 'potential_adducts' param) within a charge variant. E.g. setting this to 2 will not allow an adduct composition of '1(H+),3(Na+)' if Na+ is the least probable adduct
    type: long?
  algorithm__FeatureDeconvolution__min_rt_overlap:
    doc: Minimum overlap of the convex hull' RT intersection measured against the union from two features (if CHs are given)
    type: double?
  algorithm__FeatureDeconvolution__intensity_filter:
    doc: Enable the intensity filter, which will only allow edges between two equally charged features if the intensity of the feature with less likely adducts is smaller than that of the other feature. It is not used for features of different charge.
    type: boolean?
  algorithm__FeatureDeconvolution__negative_mode:
    doc: Enable negative ionization mode.
    type: string?
  algorithm__FeatureDeconvolution__default_map_label:
    doc: Label of map in output consensus file where all features are put by default
    type: string?
  algorithm__FeatureDeconvolution__verbose_level:
    doc: Amount of debug information given during processing.
    type: long?
outputs:
  out_cm:
    type: File
    outputBinding:
      glob: $(inputs.out_cm)
  out_fm:
    type: File
    outputBinding:
      glob: $(inputs.out_fm)
  outpairs:
    type: File
    outputBinding:
      glob: $(inputs.outpairs)
label: Decharger
doc: Decharges and merges different feature charge variants of the same peptide.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - Decharger
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json