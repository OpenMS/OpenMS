# Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: Apache-2.0
label: IDDecoyProbability
doc: "Estimates peptide probabilities using a decoy search strategy.\nWARNING: This util is deprecated."
inputs:
  in:
    doc: Identification input of combined forward decoy search (reindex with PeptideIndexer first)
    type: File?
  fwd_in:
    doc: Identification input of forward run
    type: File?
  rev_in:
    doc: Identification input of decoy run
    type: File?
  out:
    doc: Identification output with forward scores converted to probabilities
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
  decoy_algorithm__number_of_bins:
    doc: Number of bins used for the fitting, if sparse datasets are used, this number should be smaller
    type: long?
  decoy_algorithm__lower_score_better_default_value_if_zero:
    doc: This value is used if e.g. a E-value score is 0 and cannot be transformed in a real number (log of E-value)
    type: double?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - IDDecoyProbability
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json
