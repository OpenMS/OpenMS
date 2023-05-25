inputs:
  executable:
    doc: Path to the SpectraST executable to use; may be empty if the executable is globally available.
    type: File
  spectra_files:
    doc: File names(s) of spectra to be searched.
    type: File[]
  library_file:
    doc: Specify library file.
    type: File
  sequence_database_file:
    doc: The sequence database.
    type: File?
  sequence_database_type:
    doc: Specify type of sequence database
    type: string?
  search_file:
    doc: Only search a subset of the query spectra in the search file
    type: File?
  params_file:
    doc: Read search options from file. All options set in the file will be overridden by command-line options, if specified.
    type: File?
  precursor_mz_tolerance:
    doc: m/z (in Th) tolerance within which candidate entries are compared to the query. Monoisotopic mass is assumed.
    type: double?
  use_isotopically_averaged_mass:
    doc: Use isotopically averaged mass instead of monoisotopic mass
    type: boolean?
  use_all_charge_states:
    doc: Search library spectra of all charge states, i.e., ignore specified charge state (if any) of the query spectrum
    type: boolean?
  user_mod_file:
    doc: Specify name of user-defined modifications file. Default is "spectrast.usermods".
    type: File?
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
label: SpectraSTSearchAdapter
doc: Interface to the SEARCH Mode of the SpectraST executable
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - SpectraSTSearchAdapter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json