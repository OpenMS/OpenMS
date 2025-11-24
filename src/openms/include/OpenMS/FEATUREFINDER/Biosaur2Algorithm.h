// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <vector>

namespace OpenMS
{

/**
  @brief Implementation of the Biosaur2 feature detection workflow.

  The algorithm converts centroided (optionally profile) MS1 experiments into peptide features by building RT-contiguous signal traces ("hills"),
  splitting them at valleys, annotating isotopic patterns, and finally transferring the resulting characteristics into an OpenMS FeatureMap.
  Detailed descriptions of the individual processing stages as well as tunable parameters are available through @ref TOPP_Biosaur2.

  @ingroup FeatureFinder
*/
class OPENMS_DLLAPI Biosaur2Algorithm :
  public DefaultParamHandler
{
public:
  /**
    @brief Representation of a single hill (continuous m/z trace across adjacent scans).

    Each hill keeps track of spectrum/peak indices, RT/mz/intensity tuples, as well as optional mobility annotations and precomputed summary statistics.
  */
  struct Hill
  {
    std::vector<Size> scan_indices;
    std::vector<Size> peak_indices;
    std::vector<double> mz_values;
    std::vector<double> intensities;
    std::vector<double> rt_values;
    std::vector<double> drift_times;
    std::vector<double> ion_mobilities;
    double mz_median = 0.0;
    double rt_start = 0.0;
    double rt_end = 0.0;
    double rt_apex = 0.0;
    double intensity_apex = 0.0;
    double intensity_sum = 0.0;
    double drift_time_median = -1.0;
    double ion_mobility_median = -1.0;
    Size length = 0;
    Size hill_idx = 0;
  };

  /**
    @brief Candidate isotope peak that can be associated with the monoisotopic hill.

    Stores the index of the underlying hill, the ordinal isotope number (1 == first isotope), and diagnostic metrics (mass accuracy and cosine correlation).
  */
  struct IsotopeCandidate
  {
    Size hill_idx = 0;
    Size isotope_number = 0;
    double mass_diff_ppm = 0.0;
    double cos_corr = 0.0;
  };

  /**
    @brief Aggregated properties of a detected peptide feature.

    The structure retains the primary apex/shape information as well as optional FAIMS/ion-mobility annotations and the contributing isotope hills.
  */
  struct PeptideFeature
  {
    double mz = 0.0;
    double rt_start = 0.0;
    double rt_end = 0.0;
    double rt_apex = 0.0;
    double intensity_apex = 0.0;
    double intensity_sum = 0.0;
    int charge = 0;
    Size n_isotopes = 0;
    Size n_scans = 0;
    double mass_calib = 0.0;
    double drift_time = -1.0;
    double ion_mobility = -1.0;
    std::vector<IsotopeCandidate> isotopes;
    Size mono_hill_idx = 0;
  };

  Biosaur2Algorithm();

  /// @brief Set the MS data used for feature detection (copy version)
  void setMSData(const MSExperiment& ms_data);
  
  /// @brief Set the MS data used for feature detection (move version)
  void setMSData(MSExperiment&& ms_data);
  
  /// @brief Get non-const reference to MS data
  MSExperiment& getMSData();
  
  /// @brief Get const reference to MS data
  const MSExperiment& getMSData() const;

  /// Convenience overload discarding the intermediate hills/peptide-feature vectors.
  void run(FeatureMap& feature_map);

  /**
    @brief Execute the Biosaur2 workflow on the stored MS1 experiment.

    @param feature_map   Resulting FeatureMap that receives the detected features and meta information.
    @param hills         Container for all temporary hills that survived the filtering steps.
    @param peptide_features  Container storing the intermediate peptide feature representations before they are converted into FeatureMap entries.
  */
  void run(FeatureMap& feature_map,
           std::vector<Hill>& hills,
           std::vector<PeptideFeature>& peptide_features);

  /**
    @brief Export detected peptide features to a Biosaur2-compatible TSV file.

    @param features  Peptide features returned by @ref detectIsotopePatterns.
    @param filename  Destination file path.
  */
  void writeTSV(const std::vector<PeptideFeature>& features, const String& filename) const;

  /**
    @brief Export the detected hills as TSV for diagnostic purposes.

    @param hills     Hills returned by @ref detectHills / @ref processHills.
    @param filename  Destination file path.
  */
  void writeHills(const std::vector<Hill>& hills, const String& filename) const;

protected:
  void updateMembers_() override;

private:
  /// Calculate the ppm difference between two m/z values.
  double calculatePPM_(double mz1, double mz2) const;
  /// Utility median helper used for hill meta calculations.
  double calculateMedian_(const std::vector<double>& values) const;
  /// Apply a simple mean filter with the specified window size.
  std::vector<double> meanFilter_(const std::vector<double>& data, Size window) const;
  /// Estimate shift/sigma from a collection of mass errors (used for hill/isotope calibration).
  std::pair<double, double> calibrateMass_(const std::vector<double>& mass_errors, double bin_width = 0.05) const;
  /// Perform TOF-specific intensity thresholding prior to hill detection.
  void processTOF_(MSExperiment& exp) const;
  /// Centroid profile spectra using PeakPickerHiRes if requested via parameters.
  void centroidProfileSpectra_(MSExperiment& exp) const;
  /// Centroid PASEF/TIMS spectra in joint m/z–ion-mobility space (2D clustering).
  void centroidPASEFData_(MSExperiment& exp, double mz_step) const;
  /// Detect hills in the input experiment and optionally collect mass differences for calibration.
  std::vector<Hill> detectHills_(const MSExperiment& exp, double htol_ppm, double min_intensity, double min_mz, double max_mz, std::vector<double>* hill_mass_diffs = nullptr) const;
  /// Filter/process hills (length checks, summary statistics).
  std::vector<Hill> processHills_(const std::vector<Hill>& hills, Size min_length) const;
  /// Split hills at valley positions based on the configured hvf parameter.
  std::vector<Hill> splitHills_(const std::vector<Hill>& hills, double hvf, Size min_length) const;
  /// Evaluate whether isotope traces should be truncated at valleys.
  Size checkIsotopeValleySplit_(const std::vector<IsotopeCandidate>& isotopes, const std::vector<Hill>& hills, double ivf) const;
  /// Detect isotope patterns and build peptide features from the available hills.
  std::vector<PeptideFeature> detectIsotopePatterns_(std::vector<Hill>& hills, double itol_ppm, int min_charge, int max_charge, bool negative_mode, double ivf, int iuse, bool enable_isotope_calib) const;
  /// Transfer peptide feature structures into an OpenMS FeatureMap instance and attach convex hulls
  /// derived from the contributing hills (monoisotopic trace plus associated isotope hills).
  FeatureMap convertToFeatureMap_(const std::vector<PeptideFeature>& features,
                                  const std::vector<Hill>& hills) const;
  /// Utility to compute cosine correlation between two intensity traces (shared scan indices).
  double cosineCorrelation_(const std::vector<double>& intensities1, const std::vector<Size>& scans1,
                            const std::vector<double>& intensities2, const std::vector<Size>& scans2) const;
  /// Check if missing ion mobility array should be treated as an error for the given spectrum.
  bool shouldThrowForMissingIM_(const MSSpectrum& spectrum) const;

  /// Build FAIMS-aware processing groups (one per CV plus optional non-FAIMS group).
  std::vector<std::pair<double, MSExperiment>> buildFAIMSGroups_() const;

  /// Process a single FAIMS group: optional PASEF centroiding, hill and isotope detection.
  void processFAIMSGroup_(double faims_cv,
                          MSExperiment& group_exp,
                          double original_paseftol,
                          std::vector<Hill>& all_hills,
                          std::vector<PeptideFeature>& all_features,
                          Size& hill_idx_offset);

  MSExperiment ms_data_; ///< input LC-MS data
  
  double mini_;
  double minmz_;
  double maxmz_;
  double htol_;
  double itol_;
  double hvf_;
  double ivf_;
  Size minlh_;
  int cmin_;
  int cmax_;
  double pasefmini_;
  Size pasefminlh_;
  int threads_;
  int iuse_;
  bool negative_mode_;
  bool tof_mode_;
  bool profile_mode_;
  bool use_hill_calib_;
  bool ignore_iso_calib_;
  double paseftol_;
  double hrttol_;
};

} // namespace OpenMS
