inputs:
  type:
    doc: Isobaric Quantitation method used in the experiment.
    type: string?
  in:
    doc: "input raw/picked data file "
    type: File
  out:
    doc: output consensusXML file with quantitative information
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
  extraction__select_activation:
    doc: Operate only on MSn scans where any of its precursors features a certain activation method. Setting to "auto" uses HCD and HCID spectra. Set to empty string if you want to disable filtering.
    type: string?
  extraction__reporter_mass_shift:
    doc: Allowed shift (left to right) in Th from the expected position.
    type: double?
  extraction__min_precursor_intensity:
    doc: Minimum intensity of the precursor to be extracted. MS/MS scans having a precursor with a lower intensity will not be considered for quantitation.
    type: double?
  extraction__keep_unannotated_precursor:
    doc: Flag if precursor with missing intensity value or missing precursor spectrum should be included or not.
    type: string?
  extraction__min_reporter_intensity:
    doc: Minimum intensity of the individual reporter ions to be extracted.
    type: double?
  extraction__discard_low_intensity_quantifications:
    doc: Remove all reporter intensities if a single reporter is below the threshold given in 'min_reporter_intensity'.
    type: boolean?
  extraction__min_precursor_purity:
    doc: Minimum fraction of the total intensity in the isolation window of the precursor spectrum attributable to the selected precursor.
    type: double?
  extraction__precursor_isotope_deviation:
    doc: Maximum allowed deviation (in ppm) between theoretical and observed isotopic peaks of the precursor peak in the isolation window to be counted as part of the precursor.
    type: double?
  extraction__purity_interpolation:
    doc: If set to true the algorithm will try to compute the purity as a time weighted linear combination of the precursor scan and the following scan. If set to false, only the precursor scan will be used.
    type: string?
  itraq4plex__channel_114_description:
    doc: Description for the content of the 114 channel.
    type: string?
  itraq4plex__channel_115_description:
    doc: Description for the content of the 115 channel.
    type: string?
  itraq4plex__channel_116_description:
    doc: Description for the content of the 116 channel.
    type: string?
  itraq4plex__channel_117_description:
    doc: Description for the content of the 117 channel.
    type: string?
  itraq4plex__reference_channel:
    doc: Number of the reference channel (114-117).
    type: long?
  itraq4plex__correction_matrix:
    doc: "Correction matrix for isotope distributions (see documentation); use the following format: <-2Da>/<-1Da>/<+1Da>/<+2Da>; e.g. '0/0.3/4/0', '0.1/0.3/3/0.2'"
    type: string[]?
  itraq8plex__channel_113_description:
    doc: Description for the content of the 113 channel.
    type: string?
  itraq8plex__channel_114_description:
    doc: Description for the content of the 114 channel.
    type: string?
  itraq8plex__channel_115_description:
    doc: Description for the content of the 115 channel.
    type: string?
  itraq8plex__channel_116_description:
    doc: Description for the content of the 116 channel.
    type: string?
  itraq8plex__channel_117_description:
    doc: Description for the content of the 117 channel.
    type: string?
  itraq8plex__channel_118_description:
    doc: Description for the content of the 118 channel.
    type: string?
  itraq8plex__channel_119_description:
    doc: Description for the content of the 119 channel.
    type: string?
  itraq8plex__channel_121_description:
    doc: Description for the content of the 121 channel.
    type: string?
  itraq8plex__reference_channel:
    doc: Number of the reference channel (113-121). Please note that 120 is not valid.
    type: long?
  itraq8plex__correction_matrix:
    doc: "Correction matrix for isotope distributions (see documentation); use the following format: <-2Da>/<-1Da>/<+1Da>/<+2Da>; e.g. '0/0.3/4/0', '0.1/0.3/3/0.2'"
    type: string[]?
  quantification__isotope_correction:
    doc: Enable isotope correction (highly recommended). Note that you need to provide a correct isotope correction matrix otherwise the tool will fail or produce invalid results.
    type: string?
  quantification__normalization:
    doc: Enable normalization of channel intensities with respect to the reference channel. The normalization is done by using the Median of Ratios (every channel / Reference). Also the ratio of medians (from any channel and reference) is provided as control measure!
    type: boolean?
  tmt10plex__channel_126_description:
    doc: Description for the content of the 126 channel.
    type: string?
  tmt10plex__channel_127N_description:
    doc: Description for the content of the 127N channel.
    type: string?
  tmt10plex__channel_127C_description:
    doc: Description for the content of the 127C channel.
    type: string?
  tmt10plex__channel_128N_description:
    doc: Description for the content of the 128N channel.
    type: string?
  tmt10plex__channel_128C_description:
    doc: Description for the content of the 128C channel.
    type: string?
  tmt10plex__channel_129N_description:
    doc: Description for the content of the 129N channel.
    type: string?
  tmt10plex__channel_129C_description:
    doc: Description for the content of the 129C channel.
    type: string?
  tmt10plex__channel_130N_description:
    doc: Description for the content of the 130N channel.
    type: string?
  tmt10plex__channel_130C_description:
    doc: Description for the content of the 130C channel.
    type: string?
  tmt10plex__channel_131_description:
    doc: Description for the content of the 131 channel.
    type: string?
  tmt10plex__reference_channel:
    doc: The reference channel (126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131).
    type: string?
  tmt10plex__correction_matrix:
    doc: "Correction matrix for isotope distributions (see documentation); use the following format: <-2Da>/<-1Da>/<+1Da>/<+2Da>; e.g. '0/0.3/4/0', '0.1/0.3/3/0.2'"
    type: string[]?
  tmt11plex__channel_126_description:
    doc: Description for the content of the 126 channel.
    type: string?
  tmt11plex__channel_127N_description:
    doc: Description for the content of the 127N channel.
    type: string?
  tmt11plex__channel_127C_description:
    doc: Description for the content of the 127C channel.
    type: string?
  tmt11plex__channel_128N_description:
    doc: Description for the content of the 128N channel.
    type: string?
  tmt11plex__channel_128C_description:
    doc: Description for the content of the 128C channel.
    type: string?
  tmt11plex__channel_129N_description:
    doc: Description for the content of the 129N channel.
    type: string?
  tmt11plex__channel_129C_description:
    doc: Description for the content of the 129C channel.
    type: string?
  tmt11plex__channel_130N_description:
    doc: Description for the content of the 130N channel.
    type: string?
  tmt11plex__channel_130C_description:
    doc: Description for the content of the 130C channel.
    type: string?
  tmt11plex__channel_131N_description:
    doc: Description for the content of the 131N channel.
    type: string?
  tmt11plex__channel_131C_description:
    doc: Description for the content of the 131C channel.
    type: string?
  tmt11plex__reference_channel:
    doc: The reference channel (126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C).
    type: string?
  tmt11plex__correction_matrix:
    doc: "Correction matrix for isotope distributions (see documentation); use the following format: <-2Da>/<-1Da>/<+1Da>/<+2Da>; e.g. '0/0.3/4/0', '0.1/0.3/3/0.2'"
    type: string[]?
  tmt16plex__channel_126_description:
    doc: Description for the content of the 126 channel.
    type: string?
  tmt16plex__channel_127N_description:
    doc: Description for the content of the 127N channel.
    type: string?
  tmt16plex__channel_127C_description:
    doc: Description for the content of the 127C channel.
    type: string?
  tmt16plex__channel_128N_description:
    doc: Description for the content of the 128N channel.
    type: string?
  tmt16plex__channel_128C_description:
    doc: Description for the content of the 128C channel.
    type: string?
  tmt16plex__channel_129N_description:
    doc: Description for the content of the 129N channel.
    type: string?
  tmt16plex__channel_129C_description:
    doc: Description for the content of the 129C channel.
    type: string?
  tmt16plex__channel_130N_description:
    doc: Description for the content of the 130N channel.
    type: string?
  tmt16plex__channel_130C_description:
    doc: Description for the content of the 130C channel.
    type: string?
  tmt16plex__channel_131N_description:
    doc: Description for the content of the 131N channel.
    type: string?
  tmt16plex__channel_131C_description:
    doc: Description for the content of the 131C channel.
    type: string?
  tmt16plex__channel_132N_description:
    doc: Description for the content of the 132N channel.
    type: string?
  tmt16plex__channel_132C_description:
    doc: Description for the content of the 132C channel.
    type: string?
  tmt16plex__channel_133N_description:
    doc: Description for the content of the 133N channel.
    type: string?
  tmt16plex__channel_133C_description:
    doc: Description for the content of the 133C channel.
    type: string?
  tmt16plex__channel_134N_description:
    doc: Description for the content of the 134N channel.
    type: string?
  tmt16plex__reference_channel:
    doc: The reference channel (126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C, 132N, 132C, 133N, 133C, 134N).
    type: string?
  tmt16plex__correction_matrix:
    doc: "Correction matrix for isotope distributions in percent from the Thermo data sheet (see documentation); Please provide 16 entries (rows), separated by comma, where each entry contains 8 values in the following format: <-2C13>/<-N15-C13>/<-C13>/<-N15>/<+N15>/<+C13>/<+N15+C13>/<+2C13> e.g. one row may look like this: 'NA/0.00  /  0.82/0.65  /  NA/8.13  /  NA/0.26'. You may use whitespaces at your leisure to ease reading."
    type: string[]?
  tmt18plex__channel_126_description:
    doc: Description for the content of the 126 channel.
    type: string?
  tmt18plex__channel_127N_description:
    doc: Description for the content of the 127N channel.
    type: string?
  tmt18plex__channel_127C_description:
    doc: Description for the content of the 127C channel.
    type: string?
  tmt18plex__channel_128N_description:
    doc: Description for the content of the 128N channel.
    type: string?
  tmt18plex__channel_128C_description:
    doc: Description for the content of the 128C channel.
    type: string?
  tmt18plex__channel_129N_description:
    doc: Description for the content of the 129N channel.
    type: string?
  tmt18plex__channel_129C_description:
    doc: Description for the content of the 129C channel.
    type: string?
  tmt18plex__channel_130N_description:
    doc: Description for the content of the 130N channel.
    type: string?
  tmt18plex__channel_130C_description:
    doc: Description for the content of the 130C channel.
    type: string?
  tmt18plex__channel_131N_description:
    doc: Description for the content of the 131N channel.
    type: string?
  tmt18plex__channel_131C_description:
    doc: Description for the content of the 131C channel.
    type: string?
  tmt18plex__channel_132N_description:
    doc: Description for the content of the 132N channel.
    type: string?
  tmt18plex__channel_132C_description:
    doc: Description for the content of the 132C channel.
    type: string?
  tmt18plex__channel_133N_description:
    doc: Description for the content of the 133N channel.
    type: string?
  tmt18plex__channel_133C_description:
    doc: Description for the content of the 133C channel.
    type: string?
  tmt18plex__channel_134N_description:
    doc: Description for the content of the 134N channel.
    type: string?
  tmt18plex__channel_134C_description:
    doc: Description for the content of the 134C channel.
    type: string?
  tmt18plex__channel_135N_description:
    doc: Description for the content of the 135N channel.
    type: string?
  tmt18plex__reference_channel:
    doc: The reference channel (126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C, 132N, 132C, 133N, 133C, 134N, 134C, 135N).
    type: string?
  tmt18plex__correction_matrix:
    doc: "Correction matrix for isotope distributions in percent from the Thermo data sheet (see documentation); Please provide 18 entries (rows), separated by comma, where each entry contains 8 values in the following format: <-2C13>/<-N15-C13>/<-C13>/<-N15>/<+N15>/<+C13>/<+N15+C13>/<+2C13> e.g. one row may look like this: 'NA/0.00  /  0.82/0.65  /  NA/8.13  /  NA/0.26'. You may use whitespaces at your leisure to ease reading."
    type: string[]?
  tmt6plex__channel_126_description:
    doc: Description for the content of the 126 channel.
    type: string?
  tmt6plex__channel_127_description:
    doc: Description for the content of the 127 channel.
    type: string?
  tmt6plex__channel_128_description:
    doc: Description for the content of the 128 channel.
    type: string?
  tmt6plex__channel_129_description:
    doc: Description for the content of the 129 channel.
    type: string?
  tmt6plex__channel_130_description:
    doc: Description for the content of the 130 channel.
    type: string?
  tmt6plex__channel_131_description:
    doc: Description for the content of the 131 channel.
    type: string?
  tmt6plex__reference_channel:
    doc: Number of the reference channel (126-131).
    type: long?
  tmt6plex__correction_matrix:
    doc: "Correction matrix for isotope distributions (see documentation); use the following format: <-2Da>/<-1Da>/<+1Da>/<+2Da>; e.g. '0/0.3/4/0', '0.1/0.3/3/0.2'"
    type: string[]?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: IsobaricAnalyzer
doc: Calculates isobaric quantitative values for peptides
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - IsobaricAnalyzer
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json