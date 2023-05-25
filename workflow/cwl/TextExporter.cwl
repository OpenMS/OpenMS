inputs:
  in:
    doc: "Input file "
    type: File
  out:
    doc: Output file.
    type: string
  out_type:
    doc: "Output file type -- default: determined from file extension, ambiguous file extensions are interpreted as tsv"
    type: string?
  replacement:
    doc: Used to replace occurrences of the separator in strings before writing, if 'quoting' is 'none'
    type: string?
  quoting:
    doc: "Method for quoting of strings: 'none' for no quoting, 'double' for quoting with doubling of embedded quotes,\n'escape' for quoting with backslash-escaping of embedded quotes"
    type: string?
  no_ids:
    doc: Suppresses output of identification data.
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
  feature__minimal:
    doc: "Set this flag to write only three attributes: RT, m/z, and intensity."
    type: boolean?
  feature__add_metavalues:
    doc: Add columns for meta values which occur with a certain frequency (0-100%). Set to -1 to omit meta values (default).
    type: long?
  id__proteins_only:
    doc: Set this flag if you want only protein information from an idXML file
    type: boolean?
  id__peptides_only:
    doc: Set this flag if you want only peptide information from an idXML file
    type: boolean?
  id__protein_groups:
    doc: Set this flag if you want to also write indist. group information from an idXML file
    type: boolean?
  id__first_dim_rt:
    doc: If this flag is set the first_dim RT of the peptide hits will also be printed (if present).
    type: boolean?
  id__add_metavalues:
    doc: Add columns for meta values of PeptideID (=spectrum) entries which occur with a certain frequency (0-100%). Set to -1 to omit meta values (default).
    type: long?
  id__add_hit_metavalues:
    doc: Add columns for meta values of PeptideHit (=PSM) entries which occur with a certain frequency (0-100%). Set to -1 to omit meta values (default).
    type: long?
  id__add_protein_hit_metavalues:
    doc: Add columns for meta values on protein level which occur with a certain frequency (0-100%). Set to -1 to omit meta values (default).
    type: long?
  consensus__centroids:
    doc: Output file for centroids of consensus features
    type: string
  consensus__elements:
    doc: Output file for elements of consensus features
    type: string
  consensus__features:
    doc: Output file for consensus features and contained elements from all maps (writes 'nan's if elements are missing)
    type: string
  consensus__sorting_method:
    doc: "Sorting options can be combined. The precedence is: sort_by_size, sort_by_maps, sorting_method"
    type: string?
  consensus__sort_by_maps:
    doc: Apply a stable sort by the covered maps, lexicographically
    type: boolean?
  consensus__sort_by_size:
    doc: Apply a stable sort by decreasing size (i.e., the number of elements)
    type: boolean?
  consensus__add_metavalues:
    doc: Add columns for ConsensusFeature meta values.
    type: boolean?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  consensus__centroids:
    type: File
    outputBinding:
      glob: $(inputs.consensus__centroids)
  consensus__elements:
    type: File
    outputBinding:
      glob: $(inputs.consensus__elements)
  consensus__features:
    type: File
    outputBinding:
      glob: $(inputs.consensus__features)
label: TextExporter
doc: Exports various XML formats to a text file.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - TextExporter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json