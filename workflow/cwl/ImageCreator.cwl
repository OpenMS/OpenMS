inputs:
  in:
    doc: "input file "
    type: File
  in_featureXML:
    doc: "input file "
    type: File?
  out:
    doc: output file
    type: string
  out_type:
    doc: The image format. Set this if you want to force a format not reflected by the 'out' filename.
    type: string?
  rt:
    doc: Retention time range to extract
    type: string?
  mz:
    doc: Mass-to-charge range to extract
    type: string?
  width:
    doc: "Number of pixels in m/z dimension.\nIf 0, one pixel per Th."
    type: long?
  height:
    doc: "Number of pixels in RT dimension.\nIf 0, one pixel per spectrum."
    type: long?
  background_color:
    doc: "Background color e.g.: \"#FF0000\" to choose red as background color"
    type: string?
  feature_color:
    doc: "Feature color e.g.: \"#00FF00\" to choose green as feature color"
    type: string?
  gradient:
    doc: "Intensity gradient that defines colors for the range between 0 and 100.\nExample: '0,#FFFFFF;50,#FF0000;100,#000000'"
    type: string?
  max_intensity:
    doc: "Maximum peak intensity used to determine range for colors.\nIf 0, this is determined from the data."
    type: double?
  log_intensity:
    doc: Apply logarithm to intensity values
    type: boolean?
  transpose:
    doc: "Flag to transpose the resampled matrix (RT vs. m/z).\nPer default, dimensions run bottom-up in RT and left-right in m/z."
    type: boolean?
  precursors:
    doc: "Mark locations of MS2 precursors.\n"
    type: boolean?
  precursor_color:
    doc: Color for precursor marks (color code or word, e.g. 'black') (requires 'precursors' flag to be active)
    type: string?
  precursor_size:
    doc: Size of the precursor marks (requires 'precursors' flag to be active)
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
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: ImageCreator
doc: Transforms an LC-MS map into an image.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - ImageCreator
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json