# Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: Apache-2.0
label: MRMMapper
doc: MRMMapper maps measured chromatograms (mzML) and the transitions used (TraML)
inputs:
  in:
    doc: Input file containing chromatograms (converted mzXML file)
    type: File
  tr:
    doc: transition file
    type: File
  out:
    doc: Output file containing mapped chromatograms
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
  algorithm__precursor_tolerance:
    doc: Precursor tolerance when mapping (in Th)
    type: double?
  algorithm__product_tolerance:
    doc: Product tolerance when mapping (in Th)
    type: double?
  algorithm__map_multiple_assays:
    doc: Allow to map multiple assays to chromatograms and duplicate these chromatograms in the output.
    type: boolean?
  algorithm__error_on_unmapped:
    doc: Treat remaining, unmapped chromatograms as an error
    type: boolean?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MRMMapper
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json
