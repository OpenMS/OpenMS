inputs:
  in:
    doc: Input files used to generate the mzTab file.
    type: File?
  out:
    doc: Output file (mzTab)
    type: string
  first_run_inference_only:
    doc: Does the first IdentificationRun in the file only represent (protein) inference results? If so, read peptide information only from second to last runs.
    type: boolean?
  export_all_psms:
    doc: Export all PSMs instead of only the best per spectrum
    type: boolean?
  opt_columns:
    doc: Add optional columns which are not part of the mzTab standard.
    type: string[]?
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
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: MzTabExporter
doc: Exports various XML formats to an mzTab file.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MzTabExporter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json