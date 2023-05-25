inputs:
  in:
    doc: Input file(s)
    type: File[]
  id_in:
    doc: Optional idXML Input file(s) in the same order as mzML files - for Maracluster Cluster annotation
    type: File[]?
  out:
    doc: Output file in idXML format
    type: string
  consensus_out:
    doc: Consensus spectra in mzML format
    type: string
  output_directory:
    doc: Output directory for MaRaCluster original consensus output
    type: string?
  pcut:
    doc: "log(p-value) cutoff, has to be < 0.0. Default: -10.0."
    type: double?
  min_cluster_size:
    doc: minimum number of spectra in a cluster for consensus spectra
    type: long?
  maracluster_executable:
    doc: The maracluster executable. Provide a full or relative path, or make sure it can be found in your PATH environment.
    type: File
  verbose:
    doc: "Set verbosity of output: 0=no processing info, 5=all."
    type: long?
  precursor_tolerance:
    doc: Precursor monoisotopic mass tolerance
    type: double?
  precursor_tolerance_units:
    doc: tolerance_mass_units 0=ppm, 1=Da
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
  consensus_out:
    type: File
    outputBinding:
      glob: $(inputs.consensus_out)
label: MaRaClusterAdapter
doc: Facilitate input to MaRaCluster and reintegrate.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MaRaClusterAdapter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json