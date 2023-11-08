# Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: Apache-2.0
label: OpenSwathFileSplitter
doc: Splits SWATH files into n files, each containing one window.
inputs:
  in:
    doc: Input file (SWATH/DIA file)
    type: File
  outputDirectory:
    doc: Output file prefix
    type: string?
  out_qc:
    doc: Optional QC meta data (charge distribution in MS1). Only works with mzML input files.
    type: string?
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
  outputDirectory:
    type: File?
    outputBinding:
      glob: $(inputs.outputDirectory)*
  out_qc:
    type: File?
    outputBinding:
      glob: $(inputs.out_qc)
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - OpenSwathFileSplitter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json
