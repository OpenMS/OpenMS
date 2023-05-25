inputs:
  in:
    doc: Input file
    type: File
  in_type:
    doc: "Input file type -- default: determined from file extension or content"
    type: string?
  out:
    doc: Output file
    type: string
  out_type:
    doc: "Output file type -- default: determined from file extension or content"
    type: string?
  rt:
    doc: Retention time range to extract
    type: string?
  mz:
    doc: m/z range to extract (applies to ALL ms levels!)
    type: string?
  int:
    doc: Intensity range to extract
    type: string?
  sort:
    doc: Sorts the output according to RT and m/z.
    type: boolean?
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
  peak_options__sn:
    doc: Write peaks with S/N > 'sn' values only
    type: double?
  peak_options__rm_pc_charge:
    doc: Remove MS(2) spectra with these precursor charges. All spectra without precursor are kept!
    type: long[]?
  peak_options__pc_mz_range:
    doc: MSn (n>=2) precursor filtering according to their m/z value. Do not use this flag in conjunction with 'mz', unless you want to actually remove peaks in spectra (see 'mz'). RT filtering is covered by 'rt' and compatible with this flag.
    type: string?
  peak_options__pc_mz_list:
    doc: List of m/z values. If a precursor window covers ANY of these values, the corresponding MS/MS spectrum will be kept.
    type: double[]?
  peak_options__level:
    doc: MS levels to extract
    type: long[]?
  peak_options__sort_peaks:
    doc: Sorts the peaks according to m/z
    type: boolean?
  peak_options__no_chromatograms:
    doc: No conversion to space-saving real chromatograms, e.g. from SRM scans
    type: boolean?
  peak_options__remove_chromatograms:
    doc: Removes chromatograms stored in a file
    type: boolean?
  peak_options__remove_empty:
    doc: Removes spectra and chromatograms without peaks.
    type: boolean?
  peak_options__mz_precision:
    doc: Store base64 encoded m/z data using 32 or 64 bit precision
    type: string?
  peak_options__int_precision:
    doc: Store base64 encoded intensity data using 32 or 64 bit precision
    type: string?
  peak_options__indexed_file:
    doc: Whether to add an index to the file when writing
    type: string?
  peak_options__zlib_compression:
    doc: Whether to store data with zlib compression (lossless compression)
    type: boolean?
  peak_options__numpress__masstime:
    doc: "Apply MS Numpress compression algorithms in m/z or rt dimension (recommended: linear)"
    type: string?
  peak_options__numpress__lossy_mass_accuracy:
    doc: Desired (absolute) m/z accuracy for lossy compression (e.g. use 0.0001 for a mass accuracy of 0.2 ppm at 500 m/z, default uses -1.0 for maximal accuracy).
    type: double?
  peak_options__numpress__intensity:
    doc: "Apply MS Numpress compression algorithms in intensity dimension (recommended: slof or pic)"
    type: string?
  peak_options__numpress__float_da:
    doc: "Apply MS Numpress compression algorithms for the float data arrays (recommended: slof or pic)"
    type: string?
  spectra__remove_zoom:
    doc: Remove zoom (enhanced resolution) scans
    type: boolean?
  spectra__remove_mode:
    doc: Remove scans by scan mode
    type: string?
  spectra__remove_activation:
    doc: Remove MSn scans where any of its precursors features a certain activation method
    type: string?
  spectra__remove_collision_energy:
    doc: Remove MSn scans with a collision energy in the given interval
    type: string?
  spectra__remove_isolation_window_width:
    doc: Remove MSn scans whose isolation window width is in the given interval
    type: string?
  spectra__select_zoom:
    doc: Select zoom (enhanced resolution) scans
    type: boolean?
  spectra__select_mode:
    doc: "Selects scans by scan mode\n"
    type: string?
  spectra__select_activation:
    doc: Retain MSn scans where any of its precursors features a certain activation method
    type: string?
  spectra__select_collision_energy:
    doc: Select MSn scans with a collision energy in the given interval
    type: string?
  spectra__select_isolation_window_width:
    doc: Select MSn scans whose isolation window width is in the given interval
    type: string?
  spectra__select_polarity:
    doc: Retain MSn scans with a certain scan polarity
    type: string?
  spectra__replace_pc_charge:
    doc: Replaces in_charge with out_charge in all precursors.
    type: string?
  spectra__blackorwhitelist__file:
    doc: "Input file containing MS2 spectra that should be retained or removed from the mzML file!\nMatching tolerances are taken from 'spectra:blackorwhitelist:similarity_threshold|rt|mz' options.\n"
    type: File?
  spectra__blackorwhitelist__similarity_threshold:
    doc: Similarity threshold when matching MS2 spectra. (-1 = disabled).
    type: double?
  spectra__blackorwhitelist__rt:
    doc: Retention tolerance [s] when matching precursor positions. (-1 = disabled)
    type: double?
  spectra__blackorwhitelist__mz:
    doc: m/z tolerance [Th] when matching precursor positions. (-1 = disabled)
    type: double?
  spectra__blackorwhitelist__use_ppm_tolerance:
    doc: If ppm tolerance should be used. Otherwise Da are used.
    type: string?
  spectra__blackorwhitelist__blacklist:
    doc: "True: remove matched MS2. False: retain matched MS2 spectra. Other levels are kept"
    type: string?
  feature__q:
    doc: Overall quality range to extract [0:1]
    type: string?
  consensus__map:
    doc: Non-empty list of maps to be extracted from a consensus (indices are 0-based).
    type: long[]?
  consensus__map_and:
    doc: Consensus features are kept only if they contain exactly one feature from each map (as given above in 'map')
    type: boolean?
  consensus__blackorwhitelist__blacklist:
    doc: "True: remove matched MS2. False: retain matched MS2 spectra. Other levels are kept"
    type: string?
  consensus__blackorwhitelist__file:
    doc: "Input file containing consensus features whose corresponding MS2 spectra should be removed from the mzML file!\nMatching tolerances are taken from 'consensus:blackorwhitelist:rt' and 'consensus:blackorwhitelist:mz' options.\nIf consensus:blackorwhitelist:maps is specified, only these will be used.\n"
    type: File?
  consensus__blackorwhitelist__maps:
    doc: Maps used for black/white list filtering
    type: long[]?
  consensus__blackorwhitelist__rt:
    doc: Retention tolerance [s] for precursor to consensus feature position
    type: double?
  consensus__blackorwhitelist__mz:
    doc: m/z tolerance [Th] for precursor to consensus feature position
    type: double?
  consensus__blackorwhitelist__use_ppm_tolerance:
    doc: If ppm tolerance should be used. Otherwise Da are used.
    type: string?
  f_and_c__charge:
    doc: Charge range to extract
    type: string?
  f_and_c__size:
    doc: Size range to extract
    type: string?
  f_and_c__remove_meta:
    doc: Expects a 3-tuple (=3 entries in the list), i.e. <name> 'lt|eq|gt' <value>; the first is the name of meta value, followed by the comparison operator (equal, less or greater) and the value to compare to. All comparisons are done after converting the given value to the corresponding data value type of the meta value (for lists, this simply compares length, not content!)!
    type: string[]?
  f_and_c__remove_hull:
    doc: Remove hull from features.
    type: boolean?
  id__remove_clashes:
    doc: Remove features with id clashes (different sequences mapped to one feature)
    type: boolean?
  id__keep_best_score_id:
    doc: in case of multiple peptide identifications, keep only the id with best score
    type: boolean?
  id__sequences_whitelist:
    doc: Keep only features containing whitelisted substrings, e.g. features containing LYSNLVER or the modification (Oxidation). To control comparison method used for whitelisting, see 'id:sequence_comparison_method'.
    type: string[]?
  id__sequence_comparison_method:
    doc: Comparison method used to determine if a feature is whitelisted.
    type: string?
  id__accessions_whitelist:
    doc: keep only features with white listed accessions, e.g. sp|P02662|CASA1_BOVIN
    type: string[]?
  id__remove_annotated_features:
    doc: Remove features with annotations
    type: boolean?
  id__remove_unannotated_features:
    doc: Remove features without annotations
    type: boolean?
  id__remove_unassigned_ids:
    doc: Remove unassigned peptide identifications
    type: boolean?
  id__blacklist:
    doc: "Input file containing MS2 identifications whose corresponding MS2 spectra should be removed from the mzML file!\nMatching tolerances are taken from 'id:rt' and 'id:mz' options.\nThis tool will require all IDs to be matched to an MS2 spectrum, and quit with error otherwise. Use 'id:blacklist_imperfect' to allow for mismatches."
    type: File?
  id__rt:
    doc: Retention tolerance [s] for precursor to id position
    type: double?
  id__mz:
    doc: m/z tolerance [Th] for precursor to id position
    type: double?
  id__blacklist_imperfect:
    doc: Allow for mismatching precursor positions (see 'id:blacklist')
    type: boolean?
  algorithm__SignalToNoise__max_intensity:
    doc: maximal intensity considered for histogram construction. By default, it will be calculated automatically (see auto_mode). Only provide this parameter if you know what you are doing (and change 'auto_mode' to '-1')! All intensities EQUAL/ABOVE 'max_intensity' will be added to the LAST histogram bin. If you choose 'max_intensity' too small, the noise estimate might be too small as well.  If chosen too big, the bins become quite large (which you could counter by increasing 'bin_count', which increases runtime). In general, the Median-S/N estimator is more robust to a manual max_intensity than the MeanIterative-S/N.
    type: long?
  algorithm__SignalToNoise__auto_max_stdev_factor:
    doc: "parameter for 'max_intensity' estimation (if 'auto_mode' == 0): mean + 'auto_max_stdev_factor' * stdev"
    type: double?
  algorithm__SignalToNoise__auto_max_percentile:
    doc: "parameter for 'max_intensity' estimation (if 'auto_mode' == 1): auto_max_percentile th percentile"
    type: long?
  algorithm__SignalToNoise__auto_mode:
    doc: "method to use to determine maximal intensity: -1 --> use 'max_intensity'; 0 --> 'auto_max_stdev_factor' method (default); 1 --> 'auto_max_percentile' method"
    type: long?
  algorithm__SignalToNoise__win_len:
    doc: window length in Thomson
    type: double?
  algorithm__SignalToNoise__bin_count:
    doc: number of bins for intensity values
    type: long?
  algorithm__SignalToNoise__min_required_elements:
    doc: minimum number of elements required in a window (otherwise it is considered sparse)
    type: long?
  algorithm__SignalToNoise__noise_for_empty_window:
    doc: noise value used for sparse windows
    type: double?
  algorithm__SignalToNoise__write_log_messages:
    doc: Write out log messages in case of sparse windows or median in rightmost histogram bin
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: FileFilter
doc: Extracts or manipulates portions of data from peak, feature or consensus-feature files.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FileFilter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json