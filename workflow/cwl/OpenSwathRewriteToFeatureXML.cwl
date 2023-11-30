# Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: Apache-2.0
label: OpenSwathRewriteToFeatureXML
doc: Combines featureXML and mProphet tsv to FDR filtered featureXML.
inputs:
  csv:
    doc: "mProphet tsv output file: \"all_peakgroups.xls\""
    type: File?
  featureXML:
    doc: input featureXML file
    type: File
  out:
    doc: output featureXML file
    type: string
  FDR_cutoff:
    doc: FDR cutoff (e.g. to remove all features with a an m_score above 0.05 use 0.05 here)
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
  test:
    doc: Enables the test mode (needed for internal use only)
    type: boolean?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - OpenSwathRewriteToFeatureXML
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json
