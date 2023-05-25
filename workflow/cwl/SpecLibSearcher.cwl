inputs:
  in:
    doc: Input files
    type: File[]
  lib:
    doc: searchable spectral library (MSP format)
    type: File
  compare_function:
    doc: function for similarity comparison
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
  precursor__mass_tolerance:
    doc: Width of precursor mass tolerance window
    type: double?
  precursor__mass_tolerance_unit:
    doc: Unit of precursor mass tolerance.
    type: string?
  precursor__min_charge:
    doc: Minimum precursor charge to be considered.
    type: long?
  precursor__max_charge:
    doc: Maximum precursor charge to be considered.
    type: long?
  precursor__isotopes:
    doc: "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak)"
    type: long[]?
  fragment__mass_tolerance:
    doc: Fragment mass tolerance
    type: double?
  report__top_hits:
    doc: Maximum number of top scoring hits per spectrum that are reported.
    type: long?
  filter__remove_peaks_below_threshold:
    doc: All peaks of a query spectrum with intensities below <threshold> will be zeroed.
    type: double?
  filter__min_peaks:
    doc: required minimum number of peaks for a query spectrum
    type: long?
  filter__max_peaks:
    doc: Use only the top <number> of peaks.
    type: long?
  filter__cut_peaks_below:
    doc: Remove all peaks which are lower than 1/<number> of the highest peaks. Default equals all peaks which are lower than 0.001 of the maximum intensity peak
    type: long?
  modifications__fixed:
    doc: Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'
    type: string[]?
  modifications__variable:
    doc: Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'
    type: string[]?
  modifications__variable_max_per_peptide:
    doc: Maximum number of residues carrying a variable modification per candidate peptide
    type: long?
outputs:
  {}
label: SpecLibSearcher
doc: Identifies peptide MS/MS spectra by spectral matching with a searchable spectral library.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - SpecLibSearcher
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json