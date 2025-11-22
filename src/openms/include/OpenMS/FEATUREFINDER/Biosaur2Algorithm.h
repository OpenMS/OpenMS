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

class OPENMS_DLLAPI Biosaur2Algorithm :
  public DefaultParamHandler
{
public:
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

  struct IsotopeCandidate
  {
    Size hill_idx = 0;
    Size isotope_number = 0;
    double mass_diff_ppm = 0.0;
    double cos_corr = 0.0;
  };

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

  /// Convenience method discarding intermediate outputs.
  void run(const MSExperiment& input, FeatureMap& feature_map);

  void run(const MSExperiment& input,
           FeatureMap& feature_map,
           std::vector<Hill>& hills,
           std::vector<PeptideFeature>& peptide_features);

  void writeTSV(const std::vector<PeptideFeature>& features, const String& filename) const;
  void writeHills(const std::vector<Hill>& hills, const String& filename) const;

protected:
  void updateMembers_() override;

private:
  double calculatePPM(double mz1, double mz2) const;
  double calculateMedian(const std::vector<double>& values) const;
  std::vector<double> meanFilter(const std::vector<double>& data, Size window) const;
  std::pair<double, double> calibrateMass(const std::vector<double>& mass_errors, double bin_width = 0.05) const;
  void processTOF(MSExperiment& exp) const;
  void centroidProfileSpectra(MSExperiment& exp) const;
  std::vector<Hill> detectHills(const MSExperiment& exp, double htol_ppm, double min_intensity, double min_mz, double max_mz, std::vector<double>* hill_mass_diffs = nullptr) const;
  std::vector<Hill> processHills(const std::vector<Hill>& hills, Size min_length) const;
  std::vector<Hill> splitHills(const std::vector<Hill>& hills, double hvf, Size min_length) const;
  Size checkIsotopeValleySplit(const std::vector<IsotopeCandidate>& isotopes, const std::vector<Hill>& hills, double ivf) const;
  std::vector<PeptideFeature> detectIsotopePatterns(std::vector<Hill>& hills, double itol_ppm, int min_charge, int max_charge, bool negative_mode, double ivf, int iuse, bool enable_isotope_calib) const;
  FeatureMap convertToFeatureMap(const std::vector<PeptideFeature>& features) const;
  double cosineCorrelation(const std::vector<double>& intensities1, const std::vector<Size>& scans1,
                           const std::vector<double>& intensities2, const std::vector<Size>& scans2) const;

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
  int threads_;
  int iuse_;
  bool negative_mode_;
  bool tof_mode_;
  bool profile_mode_;
  bool use_hill_calib_;
  bool ignore_iso_calib_;
};

} // namespace OpenMS
