inputs:
  in1:
    doc: first input file
    type: File
  in2:
    doc: second input file
    type: File
  ratio:
    doc: acceptable relative error. Only one of 'ratio' or 'absdiff' has to be satisfied.  Use "absdiff" to deal with cases like "zero vs. epsilon".
    type: double?
  absdiff:
    doc: "acceptable absolute difference. Only one of 'ratio' or 'absdiff' has to be satisfied. "
    type: double?
  whitelist:
    doc: Lines containing one of these strings are skipped
    type: string[]?
  matched_whitelist:
    doc: Lines where one file contains one string and the other file another string are skipped. Input is given as list of colon separated tuples, e.g. String1:String2 String3:String4
    type: string[]?
  verbose:
    doc: "set verbose level:\n0 = very quiet mode (absolutely no output)\n1 = quiet mode (no output unless differences detected)\n2 = default (include summary at end)\n3 = continue after errors\n"
    type: long?
  tab_width:
    doc: tabulator width, used for calculation of column numbers
    type: long?
  first_column:
    doc: number of first column, used for calculation of column numbers
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
  {}
label: FuzzyDiff
doc: Compares two files, tolerating numeric differences.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FuzzyDiff
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json