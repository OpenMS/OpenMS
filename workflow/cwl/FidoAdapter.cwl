inputs:
  in:
    doc: "Input: identification results"
    type: File
  out:
    doc: "Output: identification results with scored/grouped proteins"
    type: string
  fido_executable:
    doc: The Fido executable. Provide a full or relative path, or make sure it can be found in your PATH environment.
    type: File
  fidocp_executable:
    doc: The FidoChooseParameters executable. Provide a full or relative path, or make sure it can be found in your PATH environment.
    type: File
  separate_runs:
    doc: Process multiple protein identification runs in the input separately, don't merge them. Merging results in loss of descriptive information of the single protein identification runs.
    type: boolean?
  keep_zero_group:
    doc: Keep the group of proteins with estimated probability of zero, which is otherwise removed (it may be very large)
    type: boolean?
  greedy_group_resolution:
    doc: Post-process Fido output with greedy resolution of shared peptides based on the protein probabilities. Also adds the resolved ambiguity groups to output.
    type: boolean?
  no_cleanup:
    doc: Omit clean-up of peptide sequences (removal of non-letter characters, replacement of I with L)
    type: boolean?
  all_PSMs:
    doc: Consider all PSMs of each peptide, instead of only the best one
    type: boolean?
  group_level:
    doc: Perform inference on protein group level (instead of individual protein level). This will lead to higher probabilities for (bigger) protein groups.
    type: boolean?
  accuracy:
    doc: Accuracy level of start parameters. There is a trade-off between accuracy and runtime.
    type: string?
  log2_states:
    doc: Binary logarithm of the max. number of connected states in a subgraph. For a value N, subgraphs that are bigger than 2^N will be split up, sacrificing accuracy for runtime. '0' uses the default (18).
    type: long?
  log2_states_precalc:
    doc: Like 'log2_states', but allows to set a separate limit for the precalculation
    type: long?
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
  prob__protein:
    doc: Protein prior probability ('gamma' parameter)
    type: double?
  prob__peptide:
    doc: Peptide emission probability ('alpha' parameter)
    type: double?
  prob__spurious:
    doc: Spurious peptide identification probability ('beta' parameter)
    type: double?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: FidoAdapter
doc: Runs the protein inference engine Fido.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FidoAdapter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json