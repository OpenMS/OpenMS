inputs:
  in:
    doc: Input protein sequences
    type: File[]
  out:
    doc: "output: simulated MS raw (profile) data"
    type: string
  out_pm:
    doc: "output: ground-truth picked (centroided) MS data"
    type: string
  out_fm:
    doc: "output: ground-truth features"
    type: string
  out_cm:
    doc: "output: ground-truth features, grouping ESI charge variants of each parent peptide"
    type: string
  out_lcm:
    doc: "output: ground-truth features, grouping labeled variants"
    type: string
  out_cntm:
    doc: "output: ground-truth features caused by contaminants"
    type: string
  out_id:
    doc: "output: ground-truth MS2 peptide identifications"
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
  algorithm__MSSim__Digestion__enzyme:
    doc: Enzyme to use for digestion (select 'no cleavage' to skip digestion)
    type: string?
  algorithm__MSSim__Digestion__model:
    doc: The cleavage model to use for digestion. 'Trained' is based on a log likelihood model (see DOI:10.1021/pr060507u).
    type: string?
  algorithm__MSSim__Digestion__min_peptide_length:
    doc: Minimum peptide length after digestion (shorter ones will be discarded)
    type: long?
  algorithm__MSSim__Digestion__model_trained__threshold:
    doc: Model threshold for calling a cleavage. Higher values increase the number of cleavages. -2 will give no cleavages, +4 almost full cleavage.
    type: double?
  algorithm__MSSim__Digestion__model_naive__missed_cleavages:
    doc: Maximum number of missed cleavages considered. All possible resulting peptides will be created.
    type: long?
  algorithm__MSSim__RT__rt_column:
    doc: Modelling of an RT or CE column
    type: string?
  algorithm__MSSim__RT__auto_scale:
    doc: Scale predicted RT's/MT's to given 'total_gradient_time'? If 'true', for CE this means that 'CE:lenght_d', 'CE:length_total', 'CE:voltage' have no influence.
    type: string?
  algorithm__MSSim__RT__total_gradient_time:
    doc: The duration [s] of the gradient.
    type: double?
  algorithm__MSSim__RT__sampling_rate:
    doc: Time interval [s] between consecutive scans
    type: double?
  algorithm__MSSim__RT__scan_window__min:
    doc: Start of RT Scan Window [s]
    type: double?
  algorithm__MSSim__RT__scan_window__max:
    doc: End of RT Scan Window [s]
    type: double?
  algorithm__MSSim__RT__variation__feature_stddev:
    doc: Standard deviation of shift in retention time [s] from predicted model (applied to every single feature independently)
    type: long?
  algorithm__MSSim__RT__variation__affine_offset:
    doc: Global offset in retention time [s] from predicted model
    type: long?
  algorithm__MSSim__RT__variation__affine_scale:
    doc: Global scaling in retention time from predicted model
    type: long?
  algorithm__MSSim__RT__column_condition__distortion:
    doc: Distortion of the elution profiles. Good presets are 0 for a perfect elution profile, 1 for a slightly distorted elution profile etc... For trapping instruments (e.g. Orbitrap) distortion should be >4.
    type: long?
  algorithm__MSSim__RT__profile_shape__width__value:
    doc: Width of the Exponential Gaussian Hybrid distribution shape of the elution profile. This does not correspond directly to the width in [s].
    type: double?
  algorithm__MSSim__RT__profile_shape__width__variance:
    doc: "Random component of the width (set to 0 to disable randomness), i.e. scale parameter for the lorentzian variation of the variance (Note: The scale parameter has to be >= 0)."
    type: double?
  algorithm__MSSim__RT__profile_shape__skewness__value:
    doc: Asymmetric component of the EGH. Higher absolute(!) values lead to more skewness (negative values cause fronting, positive values cause tailing). Tau parameter of the EGH, i.e. time constant of the exponential decay of the Exponential Gaussian Hybrid distribution shape of the elution profile.
    type: double?
  algorithm__MSSim__RT__profile_shape__skewness__variance:
    doc: "Random component of skewness (set to 0 to disable randomness), i.e. scale parameter for the lorentzian variation of the time constant (Note: The scale parameter has to be > 0)."
    type: double?
  algorithm__MSSim__RT__HPLC__model_file:
    doc: SVM model for retention time prediction
    type: string?
  algorithm__MSSim__RT__CE__pH:
    doc: pH of buffer
    type: double?
  algorithm__MSSim__RT__CE__alpha:
    doc: Exponent Alpha used to calculate mobility
    type: double?
  algorithm__MSSim__RT__CE__mu_eo:
    doc: Electroosmotic flow
    type: double?
  algorithm__MSSim__RT__CE__lenght_d:
    doc: Length of capillary [cm] from injection site to MS
    type: double?
  algorithm__MSSim__RT__CE__length_total:
    doc: Total length of capillary [cm]
    type: double?
  algorithm__MSSim__RT__CE__voltage:
    doc: Voltage applied to capillary
    type: double?
  algorithm__MSSim__Detectability__dt_simulation_on:
    doc: Modelling detectibility enabled? This can serve as a filter to remove peptides which ionize badly, thus reducing peptide count
    type: boolean?
  algorithm__MSSim__Detectability__min_detect:
    doc: Minimum peptide detectability accepted. Peptides with a lower score will be removed
    type: double?
  algorithm__MSSim__Detectability__dt_model_file:
    doc: SVM model for peptide detectability prediction
    type: string?
  algorithm__MSSim__Ionization__esi__ionized_residues:
    doc: List of residues (as three letter code) that will be considered during ES ionization. The N-term is always assumed to carry a charge. This parameter will be ignored during MALDI ionization
    type: string[]?
  algorithm__MSSim__Ionization__esi__charge_impurity:
    doc: List of charged ions that contribute to charge with weight of occurrence (their sum is scaled to 1 internally), e.g. ['H:1'] or ['H:0.7' 'Na:0.3'], ['H:4' 'Na:1'] (which internally translates to ['H:0.8' 'Na:0.2'])
    type: string[]?
  algorithm__MSSim__Ionization__esi__max_impurity_set_size:
    doc: "Maximal #combinations of charge impurities allowed (each generating one feature) per charge state. E.g. assuming charge=3 and this parameter is 2, then we could choose to allow '3H+, 2H+Na+' features (given a certain 'charge_impurity' constraints), but no '3H+, 2H+Na+, 3Na+'"
    type: long?
  algorithm__MSSim__Ionization__esi__ionization_probability:
    doc: Probability for the binomial distribution of the ESI charge states
    type: double?
  algorithm__MSSim__Ionization__maldi__ionization_probabilities:
    doc: List of probabilities for different charge states (starting at charge=1, 2, ...) during MALDI ionization (the list must sum up to 1.0)
    type: double[]?
  algorithm__MSSim__Ionization__mz__lower_measurement_limit:
    doc: Lower m/z detector limit
    type: double?
  algorithm__MSSim__Ionization__mz__upper_measurement_limit:
    doc: Upper m/z detector limit
    type: double?
  algorithm__MSSim__RawSignal__enabled:
    doc: Enable RAW signal simulation? (select 'false' if you only need feature-maps)
    type: string?
  algorithm__MSSim__RawSignal__peak_shape:
    doc: Peak Shape used around each isotope peak (be aware that the area under the curve is constant for both types, but the maximal height will differ (~ 2:3 = Lorentz:Gaussian) due to the wider base of the Lorentzian
    type: string?
  algorithm__MSSim__RawSignal__resolution__value:
    doc: Instrument resolution at 400 Th
    type: long?
  algorithm__MSSim__RawSignal__resolution__type:
    doc: How does resolution change with increasing m/z?! QTOFs usually show 'constant' behavior, FTs have linear degradation, and on Orbitraps the resolution decreases with square root of mass
    type: string?
  algorithm__MSSim__RawSignal__baseline__scaling:
    doc: Scale of baseline. Set to 0 to disable simulation of baseline
    type: double?
  algorithm__MSSim__RawSignal__baseline__shape:
    doc: The baseline is modeled by an exponential probability density function (pdf) with f(x) = shape*e^(- shape*x)
    type: double?
  algorithm__MSSim__RawSignal__mz__sampling_points:
    doc: Number of raw data points per FWHM of the peak
    type: long?
  algorithm__MSSim__RawSignal__contaminants__file:
    doc: Contaminants file with sum formula and absolute RT interval. See 'share/OpenMS/SIMULATION/contaminants.txt' for details
    type: string?
  algorithm__MSSim__RawSignal__variation__mz__error_mean:
    doc: Average systematic m/z error (in Da)
    type: double?
  algorithm__MSSim__RawSignal__variation__mz__error_stddev:
    doc: Standard deviation for m/z errors. Set to 0 to disable simulation of m/z errors
    type: double?
  algorithm__MSSim__RawSignal__variation__intensity__scale:
    doc: Constant scale factor of the feature intensity. Set to 1.0 to get the real intensity values provided in the FASTA file
    type: double?
  algorithm__MSSim__RawSignal__variation__intensity__scale_stddev:
    doc: Standard deviation of peak intensity (relative to the scaled peak height). Set to 0 to get simple rescaled intensities
    type: double?
  algorithm__MSSim__RawSignal__noise__shot__rate:
    doc: Poisson rate of shot noise per unit m/z (random peaks in m/z, where the number of peaks per unit m/z follows a Poisson distribution). Set this to 0 to disable simulation of shot noise
    type: double?
  algorithm__MSSim__RawSignal__noise__shot__intensity-mean:
    doc: Shot noise intensity mean (exponentially distributed with given mean)
    type: double?
  algorithm__MSSim__RawSignal__noise__white__mean:
    doc: Mean value of white noise (Gaussian) being added to each *measured* signal intensity
    type: double?
  algorithm__MSSim__RawSignal__noise__white__stddev:
    doc: Standard deviation of white noise being added to each *measured* signal intensity
    type: double?
  algorithm__MSSim__RawSignal__noise__detector__mean:
    doc: Mean intensity value of the detector noise (Gaussian distribution)
    type: double?
  algorithm__MSSim__RawSignal__noise__detector__stddev:
    doc: Standard deviation of the detector noise (Gaussian distribution)
    type: double?
  algorithm__MSSim__RawTandemSignal__status:
    doc: Create Tandem-MS scans?
    type: string?
  algorithm__MSSim__RawTandemSignal__tandem_mode:
    doc: "Algorithm to generate the tandem-MS spectra. 0 - fixed intensities, 1 - SVC prediction (abundant/missing), 2 - SVR prediction of peak intensity \n"
    type: long?
  algorithm__MSSim__RawTandemSignal__svm_model_set_file:
    doc: File containing the filenames of SVM Models for different charge variants
    type: string?
  algorithm__MSSim__RawTandemSignal__Precursor__ms2_spectra_per_rt_bin:
    doc: Number of allowed MS/MS spectra in a retention time bin.
    type: long?
  algorithm__MSSim__RawTandemSignal__Precursor__min_mz_peak_distance:
    doc: The minimal distance (in Th) between two peaks for concurrent selection for fragmentation. Also used to define the m/z width of an exclusion window (distance +/- from m/z of precursor). If you set this lower than the isotopic envelope of a peptide, you might get multiple fragment spectra pointing to the same precursor.
    type: double?
  algorithm__MSSim__RawTandemSignal__Precursor__mz_isolation_window:
    doc: All peaks within a mass window (in Th) of a selected peak are also selected for fragmentation.
    type: double?
  algorithm__MSSim__RawTandemSignal__Precursor__exclude_overlapping_peaks:
    doc: If true, overlapping or nearby peaks (within 'min_mz_peak_distance') are excluded for selection.
    type: boolean?
  algorithm__MSSim__RawTandemSignal__Precursor__charge_filter:
    doc: Charges considered for MS2 fragmentation.
    type: long[]?
  algorithm__MSSim__RawTandemSignal__Precursor__Exclusion__use_dynamic_exclusion:
    doc: If true dynamic exclusion is applied.
    type: boolean?
  algorithm__MSSim__RawTandemSignal__Precursor__Exclusion__exclusion_time:
    doc: The time (in seconds) a feature is excluded.
    type: double?
  algorithm__MSSim__RawTandemSignal__Precursor__ProteinBasedInclusion__max_list_size:
    doc: The maximal number of precursors in the inclusion list.
    type: long?
  algorithm__MSSim__RawTandemSignal__Precursor__ProteinBasedInclusion__rt__min_rt:
    doc: Minimal rt in seconds.
    type: double?
  algorithm__MSSim__RawTandemSignal__Precursor__ProteinBasedInclusion__rt__max_rt:
    doc: Maximal rt in seconds.
    type: double?
  algorithm__MSSim__RawTandemSignal__Precursor__ProteinBasedInclusion__rt__rt_step_size:
    doc: rt step size in seconds.
    type: double?
  algorithm__MSSim__RawTandemSignal__Precursor__ProteinBasedInclusion__rt__rt_window_size:
    doc: rt window size in seconds.
    type: long?
  algorithm__MSSim__RawTandemSignal__Precursor__ProteinBasedInclusion__thresholds__min_protein_id_probability:
    doc: Minimal protein probability for a protein to be considered identified.
    type: double?
  algorithm__MSSim__RawTandemSignal__Precursor__ProteinBasedInclusion__thresholds__min_pt_weight:
    doc: Minimal pt weight of a precursor
    type: double?
  algorithm__MSSim__RawTandemSignal__Precursor__ProteinBasedInclusion__thresholds__min_mz:
    doc: Minimal mz to be considered in protein based LP formulation.
    type: double?
  algorithm__MSSim__RawTandemSignal__Precursor__ProteinBasedInclusion__thresholds__max_mz:
    doc: Minimal mz to be considered in protein based LP formulation.
    type: double?
  algorithm__MSSim__RawTandemSignal__Precursor__ProteinBasedInclusion__thresholds__use_peptide_rule:
    doc: Use peptide rule instead of minimal protein id probability
    type: boolean?
  algorithm__MSSim__RawTandemSignal__Precursor__ProteinBasedInclusion__thresholds__min_peptide_ids:
    doc: If use_peptide_rule is true, this parameter sets the minimal number of peptide ids for a protein id
    type: long?
  algorithm__MSSim__RawTandemSignal__Precursor__ProteinBasedInclusion__thresholds__min_peptide_probability:
    doc: If use_peptide_rule is true, this parameter sets the minimal probability for a peptide to be safely identified
    type: double?
  algorithm__MSSim__RawTandemSignal__MS_E__add_single_spectra:
    doc: If true, the MS2 spectra for each peptide signal are included in the output (might be a lot). They will have a meta value 'MSE_DebugSpectrum' attached, so they can be filtered out. Native MS_E spectra will have 'MSE_Spectrum' instead.
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__isotope_model:
    doc: Model to use for isotopic peaks ('none' means no isotopic peaks are added, 'coarse' adds isotopic peaks in unit mass distance, 'fine' uses the hyperfine isotopic generator to add accurate isotopic peaks. Note that adding isotopic peaks is very slow.
    type: string?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__max_isotope:
    doc: Defines the maximal isotopic peak which is added if 'isotope_model' is 'coarse'
    type: long?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__max_isotope_probability:
    doc: Defines the maximal isotopic probability to cover if 'isotope_model' is 'fine'
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__add_metainfo:
    doc: Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__add_losses:
    doc: Adds common losses to those ion expect to have them, only water and ammonia loss is considered
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__sort_by_position:
    doc: Sort output by position
    type: string?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__add_precursor_peaks:
    doc: Adds peaks of the unfragmented precursor ion to the spectrum
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__add_all_precursor_charges:
    doc: Adds precursor peaks with all charges in the given range
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__add_abundant_immonium_ions:
    doc: Add most abundant immonium ions (for Proline, Cystein, Iso/Leucine, Histidin, Phenylalanin, Tyrosine, Tryptophan)
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__add_first_prefix_ion:
    doc: If set to true e.g. b1 ions are added
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__add_y_ions:
    doc: Add peaks of y-ions to the spectrum
    type: string?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__add_b_ions:
    doc: Add peaks of b-ions to the spectrum
    type: string?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__add_a_ions:
    doc: Add peaks of a-ions to the spectrum
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__add_c_ions:
    doc: Add peaks of c-ions to the spectrum
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__add_x_ions:
    doc: Add peaks of  x-ions to the spectrum
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__add_z_ions:
    doc: Add peaks of z-ions to the spectrum
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__y_intensity:
    doc: Intensity of the y-ions
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__b_intensity:
    doc: Intensity of the b-ions
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__a_intensity:
    doc: Intensity of the a-ions
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__c_intensity:
    doc: Intensity of the c-ions
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__x_intensity:
    doc: Intensity of the x-ions
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__z_intensity:
    doc: Intensity of the z-ions
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__relative_loss_intensity:
    doc: Intensity of loss ions, in relation to the intact ion intensity
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__precursor_intensity:
    doc: Intensity of the precursor peak
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__precursor_H2O_intensity:
    doc: Intensity of the H2O loss peak of the precursor
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__Simple__precursor_NH3_intensity:
    doc: Intensity of the NH3 loss peak of the precursor
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__add_isotopes:
    doc: If set to 1 isotope peaks of the product ion peaks are added
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__max_isotope:
    doc: Defines the maximal isotopic peak which is added, add_isotopes must be set to 1
    type: long?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__add_metainfo:
    doc: Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__add_first_prefix_ion:
    doc: If set to true e.g. b1 ions are added
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__hide_y_ions:
    doc: Add peaks of y-ions to the spectrum
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__hide_y2_ions:
    doc: Add peaks of y-ions to the spectrum
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__hide_b_ions:
    doc: Add peaks of b-ions to the spectrum
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__hide_b2_ions:
    doc: Add peaks of b-ions to the spectrum
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__hide_a_ions:
    doc: Add peaks of a-ions to the spectrum
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__hide_c_ions:
    doc: Add peaks of c-ions to the spectrum
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__hide_x_ions:
    doc: Add peaks of  x-ions to the spectrum
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__hide_z_ions:
    doc: Add peaks of z-ions to the spectrum
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__hide_losses:
    doc: Adds common losses to those ion expect to have them, only water and ammonia loss is considered
    type: boolean?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__y_intensity:
    doc: Intensity of the y-ions
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__b_intensity:
    doc: Intensity of the b-ions
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__a_intensity:
    doc: Intensity of the a-ions
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__c_intensity:
    doc: Intensity of the c-ions
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__x_intensity:
    doc: Intensity of the x-ions
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__z_intensity:
    doc: Intensity of the z-ions
    type: double?
  algorithm__MSSim__RawTandemSignal__TandemSim__SVM__relative_loss_intensity:
    doc: Intensity of loss ions, in relation to the intact ion intensity
    type: double?
  algorithm__MSSim__Global__ionization_type:
    doc: Type of Ionization (MALDI or ESI)
    type: string?
  algorithm__MSSim__Labeling__type:
    doc: Select the labeling type you want for your experiment
    type: string?
  algorithm__MSSim__Labeling__ICPL__ICPL_fixed_rtshift:
    doc: Fixed retention time shift between labeled pairs. If set to 0.0 only the retention times, computed by the RT model step are used.
    type: double?
  algorithm__MSSim__Labeling__ICPL__label_proteins:
    doc: Enables protein-labeling. (select 'false' if you only need peptide-labeling)
    type: string?
  algorithm__MSSim__Labeling__ICPL__ICPL_light_channel_label:
    doc: UniMod Id of the light channel ICPL label.
    type: string?
  algorithm__MSSim__Labeling__ICPL__ICPL_medium_channel_label:
    doc: UniMod Id of the medium channel ICPL label.
    type: string?
  algorithm__MSSim__Labeling__ICPL__ICPL_heavy_channel_label:
    doc: UniMod Id of the heavy channel ICPL label.
    type: string?
  algorithm__MSSim__Labeling__SILAC__fixed_rtshift:
    doc: Fixed retention time shift between labeled peptides. If set to 0.0 only the retention times computed by the RT model step are used.
    type: double?
  algorithm__MSSim__Labeling__SILAC__medium_channel__modification_lysine:
    doc: Modification of Lysine in the medium SILAC channel
    type: string?
  algorithm__MSSim__Labeling__SILAC__medium_channel__modification_arginine:
    doc: Modification of Arginine in the medium SILAC channel
    type: string?
  algorithm__MSSim__Labeling__SILAC__heavy_channel__modification_lysine:
    doc: Modification of Lysine in the heavy SILAC channel. If left empty, two channelSILAC is assumed.
    type: string?
  algorithm__MSSim__Labeling__SILAC__heavy_channel__modification_arginine:
    doc: Modification of Arginine in the heavy SILAC channel. If left empty, two-channel SILAC is assumed.
    type: string?
  algorithm__MSSim__Labeling__itraq__iTRAQ:
    doc: 4plex or 8plex iTRAQ?
    type: string?
  algorithm__MSSim__Labeling__itraq__reporter_mass_shift:
    doc: Allowed shift (uniformly distributed - left to right) in Da from the expected position (of e.g. 114.1, 115.1)
    type: double?
  algorithm__MSSim__Labeling__itraq__channel_active_4plex:
    doc: "Four-plex only: Each channel that was used in the experiment and its description (114-117) in format <channel>:<name>, e.g. \"114:myref\",\"115:liver\"."
    type: string[]?
  algorithm__MSSim__Labeling__itraq__channel_active_8plex:
    doc: "Eight-plex only: Each channel that was used in the experiment and its description (113-121) in format <channel>:<name>, e.g. \"113:myref\",\"115:liver\",\"118:lung\"."
    type: string[]?
  algorithm__MSSim__Labeling__itraq__isotope_correction_values_4plex:
    doc: "override default values (see Documentation); use the following format: <channel>:<-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '114:0/0.3/4/0' , '116:0.1/0.3/3/0.2' "
    type: string[]?
  algorithm__MSSim__Labeling__itraq__isotope_correction_values_8plex:
    doc: "override default values (see Documentation); use the following format: <channel>:<-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '113:0/0.3/4/0' , '116:0.1/0.3/3/0.2' "
    type: string[]?
  algorithm__MSSim__Labeling__itraq__Y_contamination:
    doc: Efficiency of labeling tyrosine ('Y') residues. 0=off, 1=full labeling
    type: double?
  algorithm__MSSim__Labeling__o18__labeling_efficiency:
    doc: Describes the distribution of the labeled peptide over the different states (unlabeled, mono- and di-labeled)
    type: double?
  algorithm__RandomNumberGenerators__biological:
    doc: Controls the 'biological' randomness of the generated data (e.g. systematic effects like deviations in RT). If set to 'random' each experiment will look different. If set to 'reproducible' each experiment will have the same outcome (given that the input data is the same)
    type: string?
  algorithm__RandomNumberGenerators__technical:
    doc: Controls the 'technical' randomness of the generated data (e.g. noise in the raw signal). If set to 'random' each experiment will look different. If set to 'reproducible' each experiment will have the same outcome (given that the input data is the same)
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_pm:
    type: File
    outputBinding:
      glob: $(inputs.out_pm)
  out_fm:
    type: File
    outputBinding:
      glob: $(inputs.out_fm)
  out_cm:
    type: File
    outputBinding:
      glob: $(inputs.out_cm)
  out_lcm:
    type: File
    outputBinding:
      glob: $(inputs.out_lcm)
  out_cntm:
    type: File
    outputBinding:
      glob: $(inputs.out_cntm)
  out_id:
    type: File
    outputBinding:
      glob: $(inputs.out_id)
label: MSSimulator
doc: A highly configurable simulator for mass spectrometry experiments.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MSSimulator
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json