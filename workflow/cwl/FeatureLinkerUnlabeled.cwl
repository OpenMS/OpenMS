# Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: Apache-2.0
label: FeatureLinkerUnlabeled
doc: Groups corresponding features from multiple maps.
inputs:
  in:
    doc: input files separated by blanks
    type: File[]
  out:
    doc: Output file
    type: string
  design:
    doc: input file containing the experimental design
    type: File?
  keep_subelements:
    doc: "For consensusXML input only: If set, the sub-features of the inputs are transferred to the output."
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
  test:
    doc: Enables the test mode (needed for internal use only)
    type: boolean?
  algorithm__second_nearest_gap:
    doc: Only link features whose distance to the second nearest neighbors (for both sides) is larger by 'second_nearest_gap' than the distance between the matched pair itself.
    type: double?
  algorithm__use_identifications:
    doc: Never link features that are annotated with different peptides (features without ID's always match; only the best hit per peptide identification is considered).
    type: boolean?
  algorithm__ignore_charge:
    doc: "false [default]: pairing requires equal charge state (or at least one unknown charge '0'); true: Pairing irrespective of charge state"
    type: boolean?
  algorithm__ignore_adduct:
    doc: "true [default]: pairing requires equal adducts (or at least one without adduct annotation); true: Pairing irrespective of adducts"
    type: string?
  algorithm__distance_RT__max_difference:
    doc: Never pair features with a larger RT distance (in seconds).
    type: double?
  algorithm__distance_RT__exponent:
    doc: Normalized RT differences ([0-1], relative to 'max_difference') are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)
    type: double?
  algorithm__distance_RT__weight:
    doc: Final RT distances are weighted by this factor
    type: double?
  algorithm__distance_MZ__max_difference:
    doc: Never pair features with larger m/z distance (unit defined by 'unit')
    type: double?
  algorithm__distance_MZ__unit:
    doc: Unit of the 'max_difference' parameter
    type: string?
  algorithm__distance_MZ__exponent:
    doc: Normalized ([0-1], relative to 'max_difference') m/z differences are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)
    type: double?
  algorithm__distance_MZ__weight:
    doc: Final m/z distances are weighted by this factor
    type: double?
  algorithm__distance_intensity__exponent:
    doc: Differences in relative intensity ([0-1]) are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)
    type: double?
  algorithm__distance_intensity__weight:
    doc: Final intensity distances are weighted by this factor
    type: double?
  algorithm__distance_intensity__log_transform:
    doc: Log-transform intensities? If disabled, d = |int_f2 - int_f1| / int_max. If enabled, d = |log(int_f2 + 1) - log(int_f1 + 1)| / log(int_max + 1))
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FeatureLinkerUnlabeled
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json
