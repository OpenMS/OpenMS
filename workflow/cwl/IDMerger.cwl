inputs:
  in:
    doc: Input files separated by blanks (all must have the same type)
    type: File[]
  out:
    doc: Output file (must have the same type as the input files)
    type: string
  out_type:
    doc: "Output file type (default: determined from file extension)"
    type: string?
  add_to:
    doc: Optional input file. IDs from 'in' are added to this file, but only if the (modified) peptide sequences are not present yet (considering only best hits per spectrum).
    type: File?
  annotate_file_origin:
    doc: "Store the original filename in each protein/peptide identification (meta value: 'file_origin') - idXML input/output only"
    type: string?
  pepxml_protxml:
    doc: "Merge idXML files derived from a pepXML and corresponding protXML file.\nExactly two input files are expected in this case. Not compatible with 'add_to'."
    type: boolean?
  merge_proteins_add_PSMs:
    doc: Merge all identified proteins by accession into one protein identification run but keep all the PSMs with updated links to potential new protein ID#s. Not compatible with 'add_to'.
    type: boolean?
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
label: IDMerger
doc: Merges several protein/peptide identification files into one file.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - IDMerger
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json