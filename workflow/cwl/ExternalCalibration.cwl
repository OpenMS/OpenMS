# Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: Apache-2.0
label: ExternalCalibration
doc: Applies an external mass recalibration.
inputs:
  in:
    doc: Input peak file
    type: File
  out:
    doc: "Output file "
    type: string
  offset:
    doc: Mass offset in ppm
    type: double?
  slope:
    doc: Slope (dependent on m/z)
    type: double?
  power:
    doc: Power (dependent on m/z)
    type: double?
  ms_level:
    doc: Target MS levels to apply the transformation onto. Scans with other levels remain unchanged.
    type: long[]?
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
  - ExternalCalibration
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json
