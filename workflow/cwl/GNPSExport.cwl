inputs:
  in_cm:
    doc: Input consensusXML file containing only consensusElements with "peptide" annotations.
    type: File
  in_mzml:
    doc: "Original mzml files containing the ms2 spectra (aka peptide annotation). \nMust be in order that the consensusXML file maps the original mzML files."
    type: File[]
  out:
    doc: Output MGF file.
    type: string
  out_quantification:
    doc: Output feature quantification table.
    type: string
  out_pairs:
    doc: Output supplementary pairs table for IIMN.
    type: string
  out_meta_values:
    doc: Output meta value file.
    type: string
  output_type:
    doc: specificity of mgf output information
    type: string?
  peptide_cutoff:
    doc: Number of most intense peptides to consider per consensus element; '-1' to consider all identifications.
    type: long?
  ms2_bin_size:
    doc: Bin size (Da) for fragment ions when merging ms2 scans.
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
  merged_spectra__cos_similarity:
    doc: Cosine similarity threshold for merged_spectra output.
    type: double?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_quantification:
    type: File
    outputBinding:
      glob: $(inputs.out_quantification)
  out_pairs:
    type: File
    outputBinding:
      glob: $(inputs.out_pairs)
  out_meta_values:
    type: File
    outputBinding:
      glob: $(inputs.out_meta_values)
label: GNPSExport
doc: "Tool to export representative consensus MS/MS scan per consensusElement into a .MGF file format.\nSee the documentation on https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking-with-openms"
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - GNPSExport
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json