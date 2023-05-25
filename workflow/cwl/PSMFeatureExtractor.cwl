inputs:
  in:
    doc: Input file(s)
    type: File[]
  out:
    doc: Output file in mzid or idXML format
    type: string
  out_type:
    doc: "Output file type -- default: determined from file extension or content."
    type: string?
  extra:
    doc: List of the MetaData parameters to be included in a feature set for precolator.
    type: string[]?
  multiple_search_engines:
    doc: Combine PSMs from different search engines by merging on scan level.
    type: boolean?
  skip_db_check:
    doc: Manual override to skip the check if same settings for multiple search engines were applied. Only valid together with -multiple_search_engines flag.
    type: boolean?
  concat:
    doc: "Naive merging of PSMs from different search engines: concatenate multiple search results instead of merging on scan level. Only valid together with -multiple_search_engines flag."
    type: boolean?
  impute:
    doc: Will instead of discarding all PSM not unanimously detected by all SE, impute missing values by their respective scores min/max observed. Only valid together with -multiple_search_engines flag.
    type: boolean?
  limit_imputation:
    doc: Will impute missing scores with the worst numerical limit (instead of min/max observed) of the respective score. Only valid together with -multiple_search_engines flag.
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
label: PSMFeatureExtractor
doc: Computes extra features for each input PSM.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - PSMFeatureExtractor
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json