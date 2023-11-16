# Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: Apache-2.0
label: OpenMSDatabasesInfo
doc: Prints the content of OpenMS' enzyme and modification databases to TSV
inputs:
  enzymes_out:
    doc: Currently supported enzymes as TSV
    type: string
  mods_out:
    doc: Currently supported modifications as TSV
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
  test:
    doc: Enables the test mode (needed for internal use only)
    type: boolean?
outputs:
  enzymes_out:
    type: File
    outputBinding:
      glob: $(inputs.enzymes_out)
  mods_out:
    type: File
    outputBinding:
      glob: $(inputs.mods_out)
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - OpenMSDatabasesInfo
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json
