inputs:
  in:
    doc: Input file (centroided data)
    type: File
  out:
    doc: Output file
    type: string
  out_csv:
    doc: "Optional CSV output file for results on 'nearest_peak' or 'highest_intensity_peak' algorithm (see corresponding subsection) containing columns: RT, uncorrectedMZ, correctedMZ, deltaMZ."
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
  feature__in:
    doc: Features used to correct precursor masses.
    type: File?
  feature__mz_tolerance:
    doc: The precursor mass tolerance. Used to determine matching to feature mass traces.
    type: double?
  feature__mz_tolerance_unit:
    doc: Unit of precursor mass tolerance
    type: string?
  feature__rt_tolerance:
    doc: Additional retention time tolerance added to feature boundaries.
    type: double?
  feature__max_trace:
    doc: Maximum isotopic trace considered in matching a precursor to a feature.
    type: long?
  feature__believe_charge:
    doc: Assume precursor charge to be correct.
    type: boolean?
  feature__keep_original:
    doc: Make a copy of the precursor and MS2 (true) or discard the original (false).
    type: boolean?
  feature__assign_all_matching:
    doc: Correct a precursor using all matching features (true) or only the nearest (false). Only evaluated if copies are created (feature:keep_original).
    type: boolean?
  nearest_peak__mz_tolerance:
    doc: The precursor mass tolerance to find the closest MS1 peak. (Disable method by setting value to 0.0)
    type: double?
  nearest_peak__mz_tolerance_unit:
    doc: Unit of precursor mass tolerance
    type: string?
  highest_intensity_peak__mz_tolerance:
    doc: The precursor mass tolerance to find the highest intensity MS1 peak. Suggested value 1/max. expected charge. (Disable method by setting value to 0.0)
    type: double?
  highest_intensity_peak__mz_tolerance_unit:
    doc: Unit of precursor mass tolerance
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_csv:
    type: File
    outputBinding:
      glob: $(inputs.out_csv)
label: HighResPrecursorMassCorrector
doc: Corrects the precursor mass and charge determined by the instrument software.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - HighResPrecursorMassCorrector
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json