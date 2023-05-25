inputs:
  in:
    doc: Input file (data annotated with identifications)
    type: File
  out:
    doc: Output file (data with one peptide identification per feature)
    type: string
  resolve_between_features:
    doc: A map may contain multiple features with both identical (possibly modified i.e. not stripped) sequence and charge state. The feature with the 'highest intensity' is very likely the most reliable one. When switched on, the filter removes the sequence annotation from the lower intensity features, thereby resolving the multiplicity. Only the most reliable features for each (possibly modified i.e. not stripped) sequence maintain annotated with this peptide sequence.
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
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: IDConflictResolver
doc: Resolves ambiguous annotations of features with peptide identifications
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - IDConflictResolver
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json