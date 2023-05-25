inputs:
  in:
    doc: Input qcml file
    type: File?
  table:
    doc: The table containing the additional qp values in the columns. First row is considered containing the header. The target run or set names/ids are indicated by column "raw data file", so each row after the header will contain the values of qps for that run. (csv without "!)
    type: File
  mapping:
    doc: The mapping of the table header to the according qp cvs, also in csv format. The first row is considered containing the headers as in the table. The second row is considered the according qp cv accessions. (csv without "!)
    type: File
  out:
    doc: Output extended qcML file
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
label: QCImporter
doc: Imports tables with quality control parameters into qcml files.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - QCImporter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json