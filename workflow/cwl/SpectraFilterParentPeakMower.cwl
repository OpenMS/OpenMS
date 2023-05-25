inputs:
  in:
    doc: "input file "
    type: File
  out:
    doc: "output file "
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
  algorithm__window_size:
    doc: The size of the m/z window where the peaks are removed, +/- window_size.
    type: double?
  algorithm__default_charge:
    doc: If the precursor has no charge set, the default charge is assumed.
    type: long?
  algorithm__clean_all_charge_states:
    doc: Set to 1 if precursor ions of all possible charge states should be removed.
    type: long?
  algorithm__consider_NH3_loss:
    doc: Whether NH3 loss peaks from the precursor should be removed.
    type: long?
  algorithm__consider_H2O_loss:
    doc: Whether H2O loss peaks from the precursor should be removed.
    type: long?
  algorithm__reduce_by_factor:
    doc: Reduce the intensities of the precursor and related ions by a given factor (set 'set_to_zero' to 0).
    type: long?
  algorithm__factor:
    doc: Factor which is used to reduce the intensities if 'reduce_by_factor' is selected.
    type: double?
  algorithm__set_to_zero:
    doc: Reduce the intensities of the precursor and related ions to zero.
    type: long?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: SpectraFilterParentPeakMower
doc: Applies thresholdfilter to peak spectra.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - SpectraFilterParentPeakMower
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json