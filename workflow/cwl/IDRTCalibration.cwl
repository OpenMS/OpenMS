# Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: Apache-2.0
label: IDRTCalibration
doc: Calibrate RTs of peptide hits linearly to standards.
inputs:
  in:
    doc: "input file "
    type: File
  out:
    doc: "output file "
    type: string
  calibrant_1_reference:
    doc: The RT of the first calibrant in the reference file.
    type: double?
  calibrant_2_reference:
    doc: The RT of the second calibrant in the reference file.
    type: double?
  calibrant_1_input:
    doc: The RT of the first calibrant in the input file. Please note that this value needs to be set. The default value -1.0 is not allowed.
    type: double?
  calibrant_2_input:
    doc: The RT of the second calibrant in the input file. Please note that this value needs to be set. The default value -1.0 is not allowed.
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
  - IDRTCalibration
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json
