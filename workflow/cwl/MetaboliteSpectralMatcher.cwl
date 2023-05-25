inputs:
  in:
    doc: Input spectra.
    type: File
  database:
    doc: Default spectral database.
    type: File
  out:
    doc: mzTab file
    type: string
  out_spectra:
    doc: Output spectra as mzML file. Can be useful to inspect the peak map after spectra merging.
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
  algorithm__prec_mass_error_value:
    doc: Error allowed for precursor ion mass.
    type: double?
  algorithm__frag_mass_error_value:
    doc: Error allowed for product ions.
    type: double?
  algorithm__mass_error_unit:
    doc: Unit of mass error (ppm or Da)
    type: string?
  algorithm__report_mode:
    doc: "Which results shall be reported: the top-three scoring ones or the best scoring one?"
    type: string?
  algorithm__ionization_mode:
    doc: Positive or negative ionization mode?
    type: string?
  algorithm__merge_spectra:
    doc: Merge MS2 spectra with the same precursor mass.
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_spectra:
    type: File
    outputBinding:
      glob: $(inputs.out_spectra)
label: MetaboliteSpectralMatcher
doc: Perform a spectral library search.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MetaboliteSpectralMatcher
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json