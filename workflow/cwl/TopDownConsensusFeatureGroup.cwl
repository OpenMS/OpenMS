# Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: Apache-2.0
label: TopDownConsensusFeatureGroup
doc: TopDownConsensusFeatureGroup from FLASHQuant
inputs:
  in:
    doc: Input tsv files to align containing feature groups (output files of FLASHQuant)
    type: File[]
  out:
    doc: Output tsv file of consensus feature groups
    type: string
  mass_tol:
    doc: Mass tolerance (Da)
    type: double?
  mass_tol_unit:
    doc: Mass tolerance unit
    type: string?
  rt_tol:
    doc: Retention time tolerance for MedianApexRetentionTime in second
    type: long?
  quant_method:
    doc: Quantity value to use from FLASHQuant result
    type: string?
  consensus_as_input:
    doc: Set it true when input files are consensus files
    type: string?
  when_duplicate:
    doc: Method to pick a mass when multiple candidates were found in the same replicate
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
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - TopDownConsensusFeatureGroup
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json
