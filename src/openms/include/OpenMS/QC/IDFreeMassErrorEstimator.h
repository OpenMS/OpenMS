// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <deque>
#include <map>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace OpenMS
{
  /**
    @brief Identification-free precursor and fragment mass-precision estimator.

    The estimator derives measurement precision from repeated precursor observations and
    likely repeated MS2 spectra without requiring peptide identifications. Precursor evidence
    uses repeated observations with the same positive charge in a bounded RT/m/z neighbourhood.
    Fragment evidence compares strong peak centres from likely repeated spectra. High-resolution
    fragment precision is estimated from unambiguous coarse-bin peak pairs with a Gaussian plus
    uniform-background mixture in ppm, avoiding dependence on a narrow fragment-match window.
    Centroided spectra use their native peaks; profile spectra use an ephemeral three-point
    quadratic fit in log-intensity space around local maxima. Input spectra are never modified.

    The returned precision describes random measurement variation. It does not recover a
    historical database-search tolerance, fixed calibration bias, or isotope-error policy.

    The class is stateful so callers can stream spectra through @ref consumeSpectrum. The
    convenience @ref compute overload resets the estimator and processes a complete experiment.

    @ingroup QC
  */
  class OPENMS_DLLAPI IDFreeMassErrorEstimator
  {
  public:
    /// Tunable estimator thresholds. Defaults match the validated prideQC implementation.
    struct OPENMS_DLLAPI Parameters
    {
      double rt_window_seconds{120.0};
      double precursor_candidate_ppm{20.0};
      double fragment_match_da{0.2};
      double fragment_low_res_match_da{0.5};
      double fragment_very_low_res_match_da{1.0};
      Size top_peaks{50};
      Size min_matched_peaks{8};
      double min_overlap_fraction{0.25};
      Size max_candidates_per_bin{32};
      Size min_spectrum_pairs{25};
      Size min_precursor_clusters{10};
      Size min_precursor_cluster_size{3};
      Size min_fragment_pairs{200};
      Size min_tolerance_pairs{200};
      Size min_tolerance_clusters{100};
      double tolerance_sigma_multiplier{6.0};
      Size min_fragment_tolerance_pairs{1000};
      Size min_fragment_tolerance_spectra{50};
      double fragment_tolerance_sigma_multiplier{6.0};
      double fragment_high_res_max_sigma_da{0.01};
      double fragment_high_res_max_sigma_ppm{10.0};
      double fragment_high_intensity_quantile{0.75};
      double fragment_low_res_min_sigma_da{0.01};
      double fragment_low_res_min_sigma_ppm{20.0};
    };

    /// Robust summary of pairwise repeated-measurement differences.
    struct OPENMS_DLLAPI RobustError
    {
      double pairwise_median{0.0};
      double pairwise_sigma{0.0};
      double single_measurement_sigma{0.0};
      double pairwise_p95_abs{0.0};
      double robust_inlier_threshold_3sigma{0.0};
      Size robust_inlier_count{0};
      Size robust_outlier_count{0};
      double robust_inlier_fraction{0.0};
      double robust_inlier_p95_abs_centered{0.0};
    };

    enum class FragmentResolutionRegime
    {
      UNAVAILABLE,
      HIGH_RESOLUTION,
      LOW_RESOLUTION
    };

    /// Precision-derived search-tolerance suggestion for one unit system.
    struct OPENMS_DLLAPI ToleranceSuggestion
    {
      double tolerance{0.0};
      double single_measurement_sigma{0.0};
      double sigma_multiplier{0.0};
      std::string unit;
      std::string confidence;
      Size support{0};
    };

    /// Estimator support/accounting counters.
    struct OPENMS_DLLAPI Diagnostics
    {
      Size precursor_eligible_ms2{0};
      Size precursor_paired_spectra{0};
      Size precursor_clusters_used{0};
      Size fragment_eligible_ms2{0};
      Size fragment_paired_spectra{0};
      Size fragment_pairs{0};
      Size fragment_mixture_pairs{0};
      Size fragment_low_res_pairs{0};
      Size fragment_very_low_res_pairs{0};
      Size fragment_centroid_spectra{0};
      Size fragment_profile_spectra{0};
      Size fragment_profile_centroids{0};
      Size fragment_zero_delta_pairs{0};
      Size fragment_mixture_signal_pairs{0};
      double fragment_zero_delta_fraction{0.0};
      double fragment_mixture_mean_ppm{0.0};
      double fragment_mixture_reference_mz{0.0};
      double fragment_mixture_signal_fraction{0.0};
      Size fragment_mixture_iterations{0};
      bool fragment_mixture_converged{false};
      bool fragment_mixture_rejected_zero_quantization{false};
      Size fragment_high_intensity_pairs{0};
      Size fragment_high_intensity_signal_pairs{0};
      double fragment_high_intensity_threshold{0.0};
      double fragment_high_intensity_signal_fraction{0.0};
      bool fragment_high_intensity_mixture_converged{false};
      bool fragment_high_intensity_fallback_used{false};
      Size excluded_fragment_profile_or_unknown{0};
      Size excluded_fragment_unknown{0};
      Size excluded_profile_peak_pick_failure{0};
      Size excluded_low_fragment_peaks{0};
      Size excluded_missing_precursor{0};
    };

    /// Current precision estimates, optional tolerance suggestions and diagnostics.
    struct OPENMS_DLLAPI Result
    {
      std::optional<RobustError> precursor_ppm;
      std::optional<RobustError> precursor_da;
      std::optional<RobustError> fragment_ppm;
      std::optional<RobustError> fragment_da;
      std::optional<ToleranceSuggestion> precursor_tolerance_ppm;
      std::optional<ToleranceSuggestion> fragment_tolerance_ppm;
      std::optional<ToleranceSuggestion> fragment_tolerance_da;
      FragmentResolutionRegime fragment_resolution_regime{FragmentResolutionRegime::UNAVAILABLE};
      double fragment_match_window_da{0.2};
      bool fragment_window_censored{false};
      Diagnostics diagnostics;
    };

    IDFreeMassErrorEstimator();
    explicit IDFreeMassErrorEstimator(const Parameters& parameters);

    /// Replace estimator parameters and clear all accumulated evidence.
    void setParameters(const Parameters& parameters);

    const Parameters& getParameters() const;

    /// Clear all accumulated evidence while retaining the current parameters.
    void reset();

    /// Consume one spectrum. Non-MS2 spectra are ignored.
    /// Charge-less MS2 spectra with a valid precursor/isolation m/z contribute
    /// fragment evidence, but not charge-specific precursor evidence.
    void consumeSpectrum(const MSSpectrum& spectrum);

    /// Reset, consume all spectra in @p experiment, and return the resulting estimate.
    Result compute(const MSExperiment& experiment);

    /// Return the current complete estimate without mutating accumulated state.
    Result getResult() const;

    /// Return only the current precursor-ppm precision summary.
    std::optional<RobustError> getPrecursorPrecisionPPM() const;

  private:
    struct PrecursorCluster
    {
      double last_rt{0.0};
      double mean_mz{0.0};
      Size count{1};
      double last_mz{0.0};
      std::vector<double> pending_da;
      std::vector<double> pending_ppm;
      bool committed{false};
    };

    struct FragmentFingerprint
    {
      double rt{0.0};
      double precursor_mz{0.0};
      std::vector<double> mz;
      std::vector<double> relative_intensity;
    };

    struct GaussianUniformFit
    {
      double mean{0.0};
      double pairwise_sigma{0.0};
      double signal_fraction{0.0};
      Size iterations{0};
      bool converged{false};
      double log_likelihood{-std::numeric_limits<double>::infinity()};
    };

    using BinKey = std::pair<Int, Int>;

    Parameters parameters_;
    std::map<BinKey, std::vector<PrecursorCluster>> precursor_bins_;
    std::map<BinKey, std::deque<FragmentFingerprint>> fragment_bins_;
    std::vector<double> precursor_errors_ppm_;
    std::vector<double> precursor_errors_da_;
    std::vector<double> fragment_errors_da_;
    std::vector<double> fragment_errors_ppm_;
    std::vector<double> fragment_mixture_errors_ppm_;
    std::vector<double> fragment_mixture_mean_mz_;
    std::vector<double> fragment_mixture_shared_relative_intensity_;
    std::vector<double> fragment_low_res_errors_da_;
    std::vector<double> fragment_low_res_errors_ppm_;
    std::vector<double> fragment_very_low_res_errors_da_;
    std::vector<double> fragment_very_low_res_errors_ppm_;
    Diagnostics diagnostics_;

    static void validateParameters_(const Parameters& parameters);
    static std::optional<RobustError> robustError_(const std::vector<double>& values);
    static std::optional<GaussianUniformFit> gaussianUniformFit_(const std::vector<double>& values);
    static RobustError errorFromGaussianUniformFit_(
      const std::vector<double>& values,
      const GaussianUniformFit& fit);
    static std::pair<std::vector<double>, std::vector<double>> centroidTopPeaks_(const MSSpectrum& spectrum, Size limit);
    static std::pair<std::vector<double>, std::vector<double>> profilePeakCenters_(const MSSpectrum& spectrum, Size limit);
    static std::pair<std::vector<double>, std::vector<double>> matchedDeltas_(
      const std::vector<double>& left_mz,
      const std::vector<double>& right_mz,
      double tolerance_da);
    static std::pair<std::vector<double>, std::vector<double>> coarseBinnedDeltas_(
      const std::vector<double>& left_mz,
      const std::vector<double>& left_relative_intensity,
      const std::vector<double>& right_mz,
      const std::vector<double>& right_relative_intensity,
      double bin_width_da,
      std::vector<double>& shared_relative_intensity);

    void commitPrecursorCluster_(PrecursorCluster& cluster);
    void consumePrecursor_(const MSSpectrum& spectrum, double precursor_mz, Int charge);
    std::pair<std::vector<double>, std::vector<double>> fragmentPeaks_(const MSSpectrum& spectrum);
    void consumeFragments_(const MSSpectrum& spectrum, double precursor_mz, Int charge);
  };
} // namespace OpenMS
