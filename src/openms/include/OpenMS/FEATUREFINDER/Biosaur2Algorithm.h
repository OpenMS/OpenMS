// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <map>
#include <vector>

namespace OpenMS
{

/**
  @brief Implementation of the Biosaur2 feature detection workflow for LC-MS1 data.

  This algorithm is a C++ reimplementation of the Biosaur2 feature detection method, originally developed in Python.
  It processes centroided LC-MS1 data (with optional profile mode support) to detect peptide features by:

  1. Building retention time-contiguous signal traces called "hills" by linking peaks across consecutive scans
  2. Splitting hills at valley points to separate co-eluting species
  3. Detecting and annotating isotopic patterns based on expected mass differences and intensity correlations
  4. Calculating comprehensive feature properties including m/z, RT, intensity, and charge state

  The algorithm supports several advanced features:
  - FAIMS compensation voltage grouping for FAIMS-enabled instruments
  - Ion mobility-aware processing for PASEF/TIMS data
  - TOF-specific intensity filtering
  - Automatic mass calibration for improved accuracy
  - Profile mode centroiding via PeakPickerHiRes

  The implementation closely mirrors the reference Python implementation to ensure reproducible results.
  All core parameters are exposed through the parameter system and can be configured via INI files
  or programmatically. For detailed parameter descriptions and usage examples, see @ref TOPP_Biosaur2.

  Reference:
  Abdrakhimov, et al. Biosaur: An open-source Python software for liquid chromatography-mass
  spectrometry peptide feature detection with ion mobility support.
  Rapid Communications in Mass Spectrometry, 2022. https://doi.org/10.1002/rcm.9045

  @htmlinclude OpenMS_Biosaur2Algorithm.parameters

  @ingroup FeatureFinder
*/
class OPENMS_DLLAPI Biosaur2Algorithm :
  public DefaultParamHandler
{
public:
  /**
    @brief Representation of a single hill (continuous m/z trace across adjacent scans).

    A hill represents a contiguous series of peaks at similar m/z values across multiple consecutive scans.
    Hills are the fundamental building blocks for feature detection. Each hill stores the raw peak data
    (indices, m/z, intensity, RT) along with optional ion mobility information for PASEF/TIMS data.
    Summary statistics (intensity-weighted mean m/z, apex RT/intensity, etc.) are precomputed for efficient downstream processing.
  */
  struct Hill
  {
    std::vector<Size> scan_indices;        ///< Indices of spectra containing peaks of this hill
    std::vector<Size> peak_indices;        ///< Indices of peaks within each spectrum
    std::vector<double> mz_values;         ///< m/z values of peaks in this hill
    std::vector<double> intensities;       ///< Intensity values of peaks in this hill
    std::vector<double> rt_values;         ///< Retention time values corresponding to each peak
    std::vector<double> drift_times;       ///< Drift time values (TIMS data), empty if not available
    std::vector<double> ion_mobilities;    ///< Ion mobility values, empty if not available
    double mz_weighted_mean = 0.0;         ///< Intensity-weighted mean m/z (hill center)
    double rt_start = 0.0;                 ///< Retention time of first peak
    double rt_end = 0.0;                   ///< Retention time of last peak
    double rt_apex = 0.0;                  ///< Retention time at maximum intensity
    double intensity_apex = 0.0;           ///< Maximum intensity value in the hill
    double intensity_sum = 0.0;            ///< Sum of all intensities in the hill
    double drift_time_median = -1.0;       ///< Median drift time (-1 if not available)
    double ion_mobility_median = -1.0;     ///< Intensity-weighted mean ion mobility (median fallback; -1 if not available)
    Size length = 0;                       ///< Number of points/peaks in this hill
    Size hill_idx = 0;                     ///< Unique identifier for this hill
  };

  /**
    @brief Candidate isotope peak that can be associated with a monoisotopic hill.

    During isotope pattern detection, candidate hills are evaluated for their suitability as isotope peaks
    of a monoisotopic hill. This structure stores the association along with quality metrics that quantify
    how well the candidate matches the expected isotope pattern in terms of mass accuracy and intensity correlation.
  */
  struct IsotopeCandidate
  {
    Size hill_idx = 0;              ///< Index of the hill that represents this isotope peak
    Size isotope_number = 0;        ///< Ordinal isotope number (0=monoisotopic, 1=first isotope, etc.)
    double mass_diff_ppm = 0.0;     ///< Mass difference to expected isotope position in ppm
    double cos_corr = 0.0;          ///< Cosine correlation between monoisotopic and isotope intensity traces
  };

  /**
    @brief Aggregated properties of a detected peptide feature.

    A peptide feature represents a complete isotope pattern detected across multiple scans. It aggregates
    information from the monoisotopic hill and its associated isotope hills, providing comprehensive
    characterization including m/z, retention time, charge state, and intensity metrics. For ion mobility
    data, drift time and ion mobility values are also included.
  */
  struct PeptideFeature
  {
    double mz = 0.0;                              ///< Monoisotopic m/z value
    double rt_start = 0.0;                        ///< Retention time of first detection
    double rt_end = 0.0;                          ///< Retention time of last detection
    double rt_apex = 0.0;                         ///< Retention time at maximum intensity
    double intensity_apex = 0.0;                  ///< Maximum intensity across the feature
    double intensity_sum = 0.0;                   ///< Sum of intensities (based on iuse parameter)
    int charge = 0;                               ///< Detected charge state
    Size n_isotopes = 0;                          ///< Number of detected isotope peaks
    Size n_scans = 0;                             ///< Number of scans contributing to the monoisotopic hill
    double mass_calib = 0.0;                      ///< Calibrated neutral mass of the feature in Da (neutral_mass = mz * charge ± charge * proton_mass; + for negative mode, − for positive mode)
    double drift_time = -1.0;                     ///< Median drift time (-1 if not available)
    double ion_mobility = -1.0;                   ///< Intensity-weighted mean ion mobility (median fallback; -1 if not available)
    std::vector<IsotopeCandidate> isotopes;       ///< List of associated isotope peaks
    Size mono_hill_idx = 0;                       ///< Index of the monoisotopic hill
  };

  /**
    @brief Default constructor

    Initializes all algorithm parameters with their default values according to the Biosaur2 specification.
  */
  Biosaur2Algorithm();

  /**
    @brief Set the MS data used for feature detection (copy version)

    @param ms_data Input MS experiment containing centroided or profile MS1 spectra
  */
  void setMSData(const MSExperiment& ms_data);

  /**
    @brief Set the MS data used for feature detection (move version)

    @param ms_data Input MS experiment containing centroided or profile MS1 spectra (will be moved)
  */
  void setMSData(MSExperiment&& ms_data);

  /**
    @brief Get non-const reference to MS data

    @return Reference to the internal MS experiment data
  */
  MSExperiment& getMSData();

  /**
    @brief Get const reference to MS data

    @return Const reference to the internal MS experiment data
  */
  const MSExperiment& getMSData() const;

  /**
    @brief Execute the Biosaur2 workflow on the stored MS1 experiment (simplified interface)

    This is a convenience overload that discards the intermediate hills and peptide feature vectors.
    Use this when you only need the final FeatureMap output.

    @param feature_map Output FeatureMap that will receive the detected features and meta information

    @note The input MS data must be set via setMSData() before calling this method
  */
  void run(FeatureMap& feature_map);

  /**
    @brief Execute the Biosaur2 workflow on the stored MS1 experiment.

    This is the main processing method that performs the complete feature detection pipeline:
    1. Filters input data to MS1 spectra only
    2. Optionally centroids profile data
    3. Applies TOF-specific intensity filtering if enabled
    4. Groups spectra by FAIMS compensation voltage if applicable
    5. Detects hills (continuous m/z traces across scans)
    6. Splits hills at valley points
    7. Detects isotope patterns and assembles peptide features
    8. Converts features to OpenMS FeatureMap format

    @param feature_map Output FeatureMap that receives the detected features and meta information
    @param hills Output container for all detected hills that survived filtering steps. Useful for diagnostics and quality control.
    @param peptide_features Output container storing intermediate peptide feature representations before conversion to FeatureMap entries

    @note The input MS data must be set via setMSData() before calling this method
    @note All spectra with MS level != 1 will be removed from the internal MS data
    @note If profile_mode is enabled, spectra will be centroided using PeakPickerHiRes
  */
  void run(FeatureMap& feature_map,
           std::vector<Hill>& hills,
           std::vector<PeptideFeature>& peptide_features);

  /**
    @brief Export detected peptide features to a Biosaur2-compatible TSV file.

    The TSV file format matches the output of the reference Python implementation, allowing
    for easy comparison and downstream processing with Biosaur2-compatible tools.

    @param features Peptide features to export (typically obtained from run())
    @param filename Destination file path for the TSV output
  */
  void writeTSV(const std::vector<PeptideFeature>& features, const String& filename) const;

  /**
    @brief Export the detected hills as TSV for diagnostic purposes.

    This method writes detailed information about each detected hill to a TSV file,
    which is useful for debugging, quality control, and understanding the feature detection process.

    @param hills Hills to export (typically obtained from run())
    @param filename Destination file path for the TSV output
  */
  void writeHills(const std::vector<Hill>& hills, const String& filename) const;

protected:
  /// @brief Update internal member variables from parameters (called automatically when parameters change)
  void updateMembers_() override;

private:
  /// @name Internal helper methods
  //@{

  /**
    @brief Calculate the mass accuracy (ppm) between two m/z values

    @param mz1 First m/z value
    @param mz2 Second m/z value (reference)
    @return Mass difference in parts per million (ppm)
  */
  double calculatePPM_(double mz1, double mz2) const;

  /**
    @brief Calculate the median of a vector of values

    @param values Input values
    @return Median value, or 0.0 if input is empty
  */
  double calculateMedian_(const std::vector<double>& values) const;

  /**
    @brief Compute a cosine correlation between two 1D intensity vectors.

    Used internally for comparing theoretical and experimental isotope intensities.
  */
  double cosineCorrelation1D_(const std::vector<double>& v1,
                             const std::vector<double>& v2) const;

  /**
    @brief Check cosine correlation for averagine-based isotope intensities.

    Returns the best correlation and the optimal truncation position in the
    experimental vector given a minimum explained averagine fraction and
    correlation threshold.
  */
  std::pair<double, Size> checkingCosCorrelationForCarbon_(const std::vector<double>& theor_full,
                                                           const std::vector<double>& exp_full,
                                                           double thresh) const;

  /**
    @brief Compute averagine-based theoretical isotope intensities.

    Uses a C-only binomial model on a 100 Da neutral-mass grid and rescales
    the resulting probabilities to the provided monoisotopic apex intensity.
    Returns the theoretical intensities and the index of the expected
    maximum-intensity isotope.
  */
  std::pair<std::vector<double>, Size> computeAveragine_(double neutral_mass,
                                                         double apex_intensity) const;

  /**
    @brief Apply a mean filter (moving average) to smooth data

    Uses zero-padding at boundaries to match NumPy's 'same' convolution mode.

    @param data Input data to be smoothed
    @param window Half-width of the filter kernel (total kernel size = 2*window + 1)
    @return Filtered data of the same length as input
  */
  std::vector<double> meanFilter_(const std::vector<double>& data, Size window) const;

  /**
    @brief Estimate mass calibration parameters from a distribution of mass errors

    Creates a histogram of mass errors and identifies the peak to determine systematic shift and spread.

    @param mass_errors Collection of mass errors (in ppm)
    @param bin_width Width of histogram bins (default: 0.05 ppm)
    @return Pair of (shift, sigma) where shift is the systematic error and sigma is the spread
  */
  std::pair<double, double> calibrateMass_(const std::vector<double>& mass_errors, double bin_width = 0.05) const;

  /**
    @brief Determine m/z binning step for hill detection.

    Performs a pre-scan over the experiment to find the maximum m/z value
    after basic filtering and derives the m/z bin width used for fast
    hill linking.
  */
  double computeHillMzStep_(const MSExperiment& exp,
                            double htol_ppm,
                            double min_intensity,
                            double min_mz,
                            double max_mz) const;

  /**
    @brief Apply TOF-specific intensity filtering

    For TOF instruments, applies a specialized filtering step to reduce noise by considering
    local intensity distributions.

    @param exp MS experiment to be filtered (modified in place)
  */
  void processTOF_(MSExperiment& exp) const;

  /**
    @brief Centroid profile spectra using PeakPickerHiRes

    @param exp MS experiment to be centroided (modified in place)
  */
  void centroidProfileSpectra_(MSExperiment& exp) const;

  /**
    @brief Centroid PASEF/TIMS spectra in joint m/z-ion mobility space

    Performs 2D clustering of peaks in the m/z and ion mobility dimensions to reduce data complexity
    for ion mobility-enabled instruments.

    @param exp MS experiment to be centroided (modified in place)
    @param mz_step m/z binning width for clustering
  */
  void centroidPASEFData_(MSExperiment& exp, double mz_step) const;

  /**
    @brief Detect hills (continuous m/z traces) in the MS experiment

    Scans through the experiment and groups peaks with similar m/z values across consecutive scans
    into hills. Optionally collects mass differences for subsequent calibration.

    @param exp Input MS experiment
    @param htol_ppm Mass tolerance in ppm for linking peaks into hills
    @param min_intensity Minimum intensity threshold for peaks
    @param min_mz Minimum m/z value to consider
    @param max_mz Maximum m/z value to consider
    @param hill_mass_diffs Optional output container for mass differences (for calibration)
    @return Vector of detected hills
  */
  std::vector<Hill> detectHills_(const MSExperiment& exp, double htol_ppm, double min_intensity, double min_mz, double max_mz, std::vector<double>* hill_mass_diffs = nullptr) const;

  /**
    @brief Link peaks in a single scan to existing hills or start new hills.

    This method implements the core hill-linking logic for one spectrum,
    updating the running hill list and the state that is carried across
    scans (previous fast m/z dictionary, ion-mobility bins, and peak-to-hill
    assignments).
  */
  void linkScanToHills_(const MSSpectrum& spectrum,
                        Size scan_idx,
                        double htol_ppm,
                        double min_intensity,
                        double min_mz,
                        double max_mz,
                        double mz_step,
                        bool use_im_global,
                        std::vector<Hill>& hills,
                        Size& hill_idx_counter,
                        std::vector<Size>& prev_peak_to_hill,
                        const MSSpectrum*& prev_spectrum_ptr,
                        std::map<int, std::vector<int>>& prev_fast_dict,
                        std::vector<int>& prev_im_bins,
                        std::vector<double>* hill_mass_diffs) const;

  /**
    @brief Filter and process hills by applying length constraints and computing summary statistics

    @param hills Input hills
    @param min_length Minimum number of scans required for a hill to be retained
    @return Filtered hills with computed statistics (median m/z, apex RT/intensity, etc.)
  */
  std::vector<Hill> processHills_(const std::vector<Hill>& hills, Size min_length) const;

  /**
    @brief Split hills at valley positions to separate co-eluting species

    Identifies local minima in the intensity profile and splits hills where the valley factor
    criterion is satisfied, effectively separating features that were initially grouped together.

    @param hills Input hills to be split
    @param hvf Hill valley factor threshold (ratio of valley to neighboring peaks)
    @param min_length Minimum length required for split segments to be retained
    @return Hills after splitting, with valley-separated segments as independent hills
  */
  std::vector<Hill> splitHills_(const std::vector<Hill>& hills, double hvf, Size min_length) const;

  /**
    @brief Evaluate whether isotope pattern should be truncated at valley positions

    Checks if there are significant valleys within the isotope traces that suggest the pattern
    should be cut short.

    @param isotopes Isotope candidates to evaluate
    @param hills All available hills
    @param ivf Isotope valley factor threshold
    @return Recommended number of isotopes to retain (may be shorter than input)
  */
  Size checkIsotopeValleySplit_(const std::vector<IsotopeCandidate>& isotopes, const std::vector<Hill>& hills, double ivf) const;

  /**
    @brief Detect isotope patterns and assemble peptide features

    For each candidate monoisotopic hill, searches for matching isotope peaks based on expected
    mass differences and charge states. Evaluates isotope candidates using cosine correlation
    and mass accuracy, then assembles complete peptide features.

    @param hills Input hills (will be modified to mark used hills)
    @param itol_ppm Mass tolerance in ppm for isotope matching
    @param min_charge Minimum charge state to consider
    @param max_charge Maximum charge state to consider
    @param negative_mode Whether to use negative ion mode (affects mass calculations)
    @param ivf Isotope valley factor for splitting patterns
    @param iuse Number of isotopes to use for intensity calculation (0=mono only, -1=all, etc.)
    @param enable_isotope_calib Whether to apply automatic mass calibration for isotopes
    @return Vector of detected peptide features
  */
  std::vector<PeptideFeature> detectIsotopePatterns_(std::vector<Hill>& hills, double itol_ppm, int min_charge, int max_charge, bool negative_mode, double ivf, int iuse, bool enable_isotope_calib) const;

  /**
    @brief Convert peptide features to OpenMS FeatureMap format

    Transfers all feature properties (m/z, RT, intensity, charge, etc.) to OpenMS Feature objects
    and constructs convex hulls from the contributing hills.
    The representation of convex hulls (full mass-trace hulls vs. a single RT–m/z bounding box)
    is controlled via the @em convex_hulls parameter.

    @param features Input peptide features
    @param hills All hills (needed to construct convex hulls)
    @return FeatureMap containing the converted features
  */
  FeatureMap convertToFeatureMap_(const std::vector<PeptideFeature>& features,
                                  const std::vector<Hill>& hills) const;

  /**
    @brief Compute cosine correlation between two intensity traces

    Evaluates similarity of intensity profiles between two hills by computing the cosine of the
    angle between their intensity vectors (considering only overlapping scan indices).

    @param intensities1 First intensity trace
    @param scans1 Scan indices for first trace
    @param intensities2 Second intensity trace
    @param scans2 Scan indices for second trace
    @return Cosine correlation coefficient [0, 1], where 1 indicates perfect correlation
  */
  double cosineCorrelation_(const std::vector<double>& intensities1, const std::vector<Size>& scans1,
                            const std::vector<double>& intensities2, const std::vector<Size>& scans2) const;

  /**
    @brief Check if missing ion mobility data should be treated as an error

    The Python reference implementation gracefully degrades when ion mobility arrays are missing.
    This method implements the same behavior by returning false to allow processing without IM data.

    @param spectrum Spectrum to check
    @return Currently always returns false (missing IM is not an error)
  */
  bool shouldThrowForMissingIM_(const MSSpectrum& spectrum) const;

  /**
    @brief Build FAIMS-aware processing groups

    Separates spectra by FAIMS compensation voltage (CV) to ensure features are detected
    within each CV group independently. Non-FAIMS spectra are placed in a separate group.

    @return Vector of (CV value, MS experiment) pairs
  */
  std::vector<std::pair<double, MSExperiment>> buildFAIMSGroups_() const;

  /**
    @brief Process a single FAIMS compensation voltage group

    Performs the complete feature detection pipeline (optional PASEF centroiding, hill detection,
    isotope pattern detection) for a single FAIMS CV group or non-FAIMS data.

    @param faims_cv FAIMS compensation voltage for this group (NaN for non-FAIMS)
    @param group_exp MS experiment for this group (modified in place)
    @param original_paseftol Original PASEF tolerance setting
    @param all_hills Accumulator for all detected hills
    @param all_features Accumulator for all detected features
    @param hill_idx_offset Current hill index offset for unique indexing
  */
  void processFAIMSGroup_(double faims_cv,
                          MSExperiment& group_exp,
                          double original_paseftol,
                          std::vector<Hill>& all_hills,
                          std::vector<PeptideFeature>& all_features,
                          Size& hill_idx_offset);
  //@}

  /// @name Member variables
  //@{
  MSExperiment ms_data_;      ///< Input LC-MS data

  // Algorithm parameters (cached from Param for performance)
  double mini_;               ///< Minimum intensity threshold
  double minmz_;              ///< Minimum m/z value
  double maxmz_;              ///< Maximum m/z value
  double htol_;               ///< Mass tolerance in ppm for hill detection
  double itol_;               ///< Mass tolerance in ppm for isotope pattern detection
  double hvf_;                ///< Hill valley factor for splitting hills at valleys
  double ivf_;                ///< Isotope valley factor for splitting isotope patterns
  Size minlh_;                ///< Minimum number of scans required for a hill
  int cmin_;                  ///< Minimum charge state to consider
  int cmax_;                  ///< Maximum charge state to consider
  double pasefmini_;          ///< Minimum intensity for PASEF/TIMS clusters after centroiding
  Size pasefminlh_;           ///< Minimum number of points per PASEF/TIMS cluster
  int iuse_;                  ///< Number of isotopes for intensity calculation (0=mono, -1=all, N=mono+N)
  bool negative_mode_;        ///< Whether to use negative ion mode
  bool tof_mode_;             ///< Whether to enable TOF-specific intensity filtering
  bool profile_mode_;         ///< Whether to centroid profile data using PeakPickerHiRes
  bool use_hill_calib_;       ///< Whether to use automatic hill mass tolerance calibration
  bool ignore_iso_calib_;     ///< Whether to disable automatic isotope mass calibration
  double paseftol_;           ///< Ion mobility tolerance for PASEF/TIMS data (0=disable)
  double hrttol_;             ///< Maximum RT difference between monoisotopic and isotope apex (0=disable)
  String convex_hull_mode_;   ///< Representation of feature convex hulls ("mass_traces" vs. "bounding_box")
  //@}
};

} // namespace OpenMS
