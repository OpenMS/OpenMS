inputs:
  in:
    doc: Input consensusXML with peptide intensities
    type: File
  in_design:
    doc: Experimental Design file
    type: File
  method:
    doc: Method used in the experiment(label free [LFQ], isobaric labeling [ISO]))
    type: string?
  msstats_bioreplicate:
    doc: Which column in the condition table should be used for MSstats 'BioReplicate'
    type: string?
  msstats_condition:
    doc: Which column in the condition table should be used for MSstats 'Condition'
    type: string?
  msstats_mixture:
    doc: Which column in the condition table should be used for MSstats 'Mixture'
    type: string?
  reannotate_filenames:
    doc: Overwrite MS file names in consensusXML
    type: File[]?
  labeled_reference_peptides:
    doc: If set, IsotopeLabelType is 'H', else 'L'
    type: boolean?
  retention_time_summarization_method:
    doc: How indistinguishable peptidoforms at different retention times should be treated. This is usually necessary for LFQ experiments and therefore defaults to 'max'. In case of TMT/iTRAQ, MSstatsTMT does the aggregation itself later and the parameter always resets to manual (i.e. is unused).
    type: string?
  out:
    doc: Input CSV file for MSstats.
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
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: MSstatsConverter
doc: Converter to input for MSstats
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MSstatsConverter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json