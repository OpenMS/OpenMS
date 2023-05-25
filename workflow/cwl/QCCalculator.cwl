inputs:
  in:
    doc: raw data input file (this is relevant if you want to look at MS1, MS2 and precursor peak information)
    type: File
  out:
    doc: Your QC file.
    type: string
  out_type:
    doc: "Output file type -- default: determined from file extension or content"
    type: string?
  label:
    doc: unique name for the run that can be used in a figure label
    type: string?
  name:
    doc: name of the person creating this mzQC file
    type: string?
  address:
    doc: contact address (mail/e-mail or phone)
    type: string?
  description:
    doc: description and comments about the mzQC file contents
    type: string?
  id:
    doc: Input idXML file containing the identifications. Your identifications will be exported in an easy-to-read format
    type: File?
  feature:
    doc: feature input file (this is relevant for most QC issues)
    type: File?
  consensus:
    doc: consensus input file (this is only used for charge state deconvoluted output. Use the consensusXML output form the DeCharger)
    type: File?
  remove_duplicate_features:
    doc: This flag should be set, if you work with a set of merged features.
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
label: QCCalculator
doc: Calculates basic quality parameters from MS experiments and subsequent analysis data as identification or feature detection.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - QCCalculator
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json