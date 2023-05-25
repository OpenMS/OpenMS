inputs:
  in:
    doc: Input files separated by blank
    type: File[]
  in_type:
    doc: "Input file type (default: determined from file extension or content)"
    type: string?
  out:
    doc: Output file
    type: string
  annotate_file_origin:
    doc: Store the original filename in each feature using meta value "file_origin" (for featureXML and consensusXML only).
    type: boolean?
  append_method:
    doc: "(ConsensusXML-only) Append quantitative information about features row-wise or column-wise.\n- 'append_rows' is usually used when the inputs come from the same MS run (e.g. caused by manual splitting or multiple algorithms on the same file)\n- 'append_cols' when you want to combine consensusXMLs from e.g. different fractions to be summarized in ProteinQuantifier or jointly exported with MzTabExporter."
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
  rt_concat__gap:
    doc: The amount of gap (in seconds) to insert between the RT ranges of different input files. RT concatenation is enabled if a value > 0 is set.
    type: double?
  raw__rt_auto:
    doc: Assign retention times automatically (integers starting at 1)
    type: boolean?
  raw__rt_custom:
    doc: List of custom retention times that are assigned to the files. The number of given retention times must be equal to the number of input files.
    type: double[]?
  raw__rt_filename:
    doc: Try to guess the retention time of a file based on the filename. This option is useful for merging DTA files, where filenames should contain the string 'rt' directly followed by a floating point number, e.g. 'my_spectrum_rt2795.15.dta'
    type: boolean?
  raw__ms_level:
    doc: If 1 or higher, this number is assigned to spectra as the MS level. This option is useful for DTA files which do not contain MS level information.
    type: long?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: FileMerger
doc: Merges several MS files into one file.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FileMerger
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json