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
  algorithm__MetaboliteFeatureDeconvolution__charge_min:
    doc: Minimal possible charge
    type: long?
  algorithm__MetaboliteFeatureDeconvolution__charge_max:
    doc: Maximal possible charge
    type: long?
  algorithm__MetaboliteFeatureDeconvolution__charge_span_max:
    doc: Maximal range of charges for a single analyte, i.e. observing q1=[5,6,7] implies span=3. Setting this to 1 will only find adduct variants of the same charge
    type: long?
  algorithm__MetaboliteFeatureDeconvolution__q_try:
    doc: Try different values of charge for each feature according to the above settings ('heuristic' [does not test all charges, just the likely ones] or 'all' ), or leave feature charge untouched ('feature').
    type: string?
  algorithm__MetaboliteFeatureDeconvolution__retention_max_diff:
    doc: Maximum allowed RT difference between any two features if their relation shall be determined
    type: double?
  algorithm__MetaboliteFeatureDeconvolution__retention_max_diff_local:
    doc: Maximum allowed RT difference between between two co-features, after adduct shifts have been accounted for (if you do not have any adduct shifts, this value should be equal to 'retention_max_diff', otherwise it should be smaller!)
    type: double?
  algorithm__MetaboliteFeatureDeconvolution__mass_max_diff:
    doc: Maximum allowed mass tolerance per feature. Defines a symmetric tolerance window around the feature. When looking at possible feature pairs, the allowed feature-wise errors are combined for consideration of possible adduct shifts. For ppm tolerances, each window is based on the respective observed feature mz (instead of putative experimental mzs causing the observed one)!
    type: double?
  algorithm__MetaboliteFeatureDeconvolution__unit:
    doc: Unit of the 'max_difference' parameter
    type: string?
  algorithm__MetaboliteFeatureDeconvolution__potential_adducts:
    doc: "Adducts used to explain mass differences in format: 'Elements:Charge(+/-/0):Probability[:RTShift[:Label]]', i.e. the number of '+' or '-' indicate the charge ('0' if neutral adduct), e.g. 'Ca:++:0.5' indicates +2. Probabilites have to be in (0,1]. The optional RTShift param indicates the expected RT shift caused by this adduct, e.g. '(2)H4H-4:0:1:-3' indicates a 4 deuterium label, which causes early elution by 3 seconds. As fifth parameter you can add a label for every feature with this adduct. This also determines the map number in the consensus file. Adduct element losses are written in the form 'H-2'. All provided adducts need to have the same charge sign or be neutral! Mixing of adducts with different charge directions is only allowed as neutral complexes. For example, 'H-1Na:0:0.05' can be used to model Sodium gains (with balancing deprotonation) in negative mode."
    type: string[]?
  algorithm__MetaboliteFeatureDeconvolution__max_neutrals:
    doc: Maximal number of neutral adducts(q=0) allowed. Add them in the 'potential_adducts' section!
    type: long?
  algorithm__MetaboliteFeatureDeconvolution__use_minority_bound:
    doc: Prune the considered adduct transitions by transition probabilities.
    type: string?
  algorithm__MetaboliteFeatureDeconvolution__max_minority_bound:
    doc: "Limits allowed adduct compositions and changes between compositions in the underlying graph optimization problem by introducing a probability-based threshold: the minority bound sets the maximum count of the least probable adduct (according to 'potential_adducts' param) within a charge variant with maximum charge only containing the most likely adduct otherwise. E.g., for 'charge_max' 4 and 'max_minority_bound' 2 with most probable adduct being H+ and least probable adduct being Na+, this will allow adduct compositions of '2(H+),2(Na+)' but not of '1(H+),3(Na+)'. Further, adduct compositions/changes less likely than '2(H+),2(Na+)' will be discarded as well."
    type: long?
  algorithm__MetaboliteFeatureDeconvolution__min_rt_overlap:
    doc: Minimum overlap of the convex hull' RT intersection measured against the union from two features (if CHs are given)
    type: double?
  algorithm__MetaboliteFeatureDeconvolution__intensity_filter:
    doc: Enable the intensity filter, which will only allow edges between two equally charged features if the intensity of the feature with less likely adducts is smaller than that of the other feature. It is not used for features of different charge.
    type: boolean?
  algorithm__MetaboliteFeatureDeconvolution__negative_mode:
    doc: Enable negative ionization mode.
    type: boolean?
  algorithm__MetaboliteFeatureDeconvolution__default_map_label:
    doc: Label of map in output consensus file where all features are put by default
    type: string?
  algorithm__MetaboliteFeatureDeconvolution__verbose_level:
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
label: MetaboliteAdductDecharger
doc: Decharges and merges different feature charge variants of the same metabolite.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MetaboliteAdductDecharger
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json