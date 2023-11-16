# Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: Apache-2.0
label: QCExporter
doc: Will extract several qp from several run/sets in a tabular format.
inputs:
  in:
    doc: Input qcml file
    type: File
  names:
    doc: The name of the target runs or sets to be exported from. If empty, from all will be exported.
    type: string[]?
  mapping:
    doc: The mapping of the exported table's headers to the according qp cvs. The first row is considered containing the headers as for the exported the table. The second row is considered the according qp cv accessions of the qp to be exported.
    type: File
  out_csv:
    doc: Output csv formatted quality parameter.
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
  out_csv:
    type: File
    outputBinding:
      glob: $(inputs.out_csv)
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - QCExporter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json
