inputs:
  in:
    doc: Input qcml file
    type: File?
  qp_att_acc:
    doc: Defines the qp cv accession of the qp to which the table/image is attached.
    type: string?
  cv_acc:
    doc: Defines the cv accession of the attachment.
    type: string
  run:
    doc: The file that defined the run under which the qp for the attachment is aggregated as mzML file. The file is only used to extract the run name from the file name.
    type: File?
  name:
    doc: If no file for the run was given (or if the target qp is contained in a set), at least a name of the target run/set containing the the qp for the attachment has to be given.
    type: string?
  plot:
    doc: If a plot image is to be attached to a qp, this has to be specified here.
    type: File?
  table:
    doc: If a table is to be attached to a qp, this has to be specified here.
    type: File?
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
label: QCEmbedder
doc: Attaches a table or an image to a given qc parameter.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - QCEmbedder
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json