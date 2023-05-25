inputs:
  in:
    doc: "input file "
    type: File
  out:
    doc: "output file "
    type: string
  pepnovo_executable:
    doc: The PepNovo executable. Provide a full or relative path, or make sure it can be found in your PATH environment.
    type: File
  model_directory:
    doc: Name of the directory where the model files are kept.
    type: string
  correct_pm:
    doc: Find optimal precursor mass and charge values.
    type: boolean?
  use_spectrum_charge:
    doc: Do not correct charge
    type: boolean?
  use_spectrum_mz:
    doc: Do not correct the precursor m/z value that appears in the file.
    type: boolean?
  no_quality_filter:
    doc: Do not remove low quality spectra.
    type: boolean?
  fragment_tolerance:
    doc: The fragment tolerance (between 0 and 0.75 Da. Set to -1.0 to use model's default setting)
    type: double?
  pm_tolerance:
    doc: The precursor mass tolerance (between 0 and 5.0 Da. Set to -1.0 to use model's default setting)
    type: double?
  model:
    doc: Name of the model that should be used
    type: string?
  digest:
    doc: Enzyme used for digestion (default TRYPSIN)
    type: string?
  tag_length:
    doc: Returns peptide sequence of the specified length (only lengths 3-6 are allowed)
    type: long?
  num_solutions:
    doc: Number of solutions to be computed
    type: long?
  fixed_modifications:
    doc: Fixed modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
    type: string[]?
  variable_modifications:
    doc: Variable modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
    type: string[]?
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
label: PepNovoAdapter
doc: Adapter to PepNovo supporting all PepNovo command line parameters. The results are converted from the PepNovo text outfile format into the idXML format.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - PepNovoAdapter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json