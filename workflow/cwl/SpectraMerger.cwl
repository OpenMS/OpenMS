inputs:
  in:
    doc: Input mzML file.
    type: File
  out:
    doc: Output mzML file with merged spectra.
    type: string
  merging_method:
    doc: Method of merging which should be used.
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
  algorithm__mz_binning_width:
    doc: minimum m/z distance for two data points (profile data) or peaks (centroided data) to be considered distinct. Closer data points or peaks will be merged.
    type: double?
  algorithm__mz_binning_width_unit:
    doc: Unit in which the distance between two data points or peaks is given.
    type: string?
  algorithm__sort_blocks:
    doc: Sort blocks by <?> before merging them (useful for precursor order)
    type: string?
  algorithm__average_gaussian__spectrum_type:
    doc: Spectrum type of the MS level to be averaged
    type: string?
  algorithm__average_gaussian__ms_level:
    doc: If set to be 0, each MS level will be merged from 1 to max. Otherwise, average spectra of this level. All other spectra remain unchanged.
    type: long?
  algorithm__average_gaussian__rt_FWHM:
    doc: FWHM of Gauss curve in seconds to be averaged over.
    type: double?
  algorithm__average_gaussian__cutoff:
    doc: Intensity cutoff for Gaussian. The Gaussian RT profile decreases from 1 at its apex to 0 at infinity. Spectra for which the intensity of the Gaussian drops below the cutoff do not contribute to the average.
    type: double?
  algorithm__average_gaussian__precursor_mass_tol:
    doc: PPM mass tolerance for precursor mass. If set, MSn (n>2) spectra of precursor masses within the tolerance are averaged.
    type: double?
  algorithm__average_gaussian__precursor_max_charge:
    doc: Possible maximum precursor ion charge. Effective only when average_gaussian:precursor_mass_tol option is active.
    type: long?
  algorithm__average_tophat__spectrum_type:
    doc: Spectrum type of the MS level to be averaged
    type: string?
  algorithm__average_tophat__ms_level:
    doc: If set to be 0, each MS level will be merged from 1 to max. Otherwise, average spectra of this level. All other spectra remain unchanged.
    type: long?
  algorithm__average_tophat__rt_range:
    doc: RT range to be averaged over, i.e. +/-(RT range)/2 from each spectrum.
    type: double?
  algorithm__average_tophat__rt_unit:
    doc: Unit for RT range.
    type: string?
  algorithm__block_method__ms_levels:
    doc: Merge spectra of this level. All spectra with other MS levels remain untouched.
    type: long[]?
  algorithm__block_method__rt_block_size:
    doc: Maximum number of scans to be summed up.
    type: long?
  algorithm__block_method__rt_max_length:
    doc: Maximum RT size of the block in seconds (0.0 = no size restriction).
    type: double?
  algorithm__precursor_method__mz_tolerance:
    doc: Max m/z distance of the precursor entries of two spectra to be merged in [Da].
    type: double?
  algorithm__precursor_method__mass_tolerance:
    doc: Max mass distance of the precursor entries of two spectra to be merged in [Da]. Active when set to a positive value.
    type: double?
  algorithm__precursor_method__rt_tolerance:
    doc: Max RT distance of the precursor entries of two spectra to be merged in [s].
    type: double?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: SpectraMerger
doc: Merges spectra (each MS level separately), increasing S/N ratios.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - SpectraMerger
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json