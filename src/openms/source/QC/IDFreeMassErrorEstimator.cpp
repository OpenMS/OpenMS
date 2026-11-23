// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/QC/IDFreeMassErrorEstimator.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/PeakTypeEstimator.h>
#include <OpenMS/MATH/StatisticFunctions.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace OpenMS
{
  namespace
  {
    struct PeakValue
    {
      double mz{0.0};
      double intensity{0.0};
    };

    std::vector<double> sortedFinite(std::vector<double> values)
    {
      std::erase_if(values, [](double value) { return !std::isfinite(value); });
      std::sort(values.begin(), values.end());
      return values;
    }

    double quantileSorted(const std::vector<double>& values, double q)
    {
      return Math::quantile(values.begin(), values.end(), q);
    }
  } // namespace

  IDFreeMassErrorEstimator::IDFreeMassErrorEstimator()
  {
    validateParameters_(parameters_);
  }

  IDFreeMassErrorEstimator::IDFreeMassErrorEstimator(const Parameters& parameters) : parameters_(parameters)
  {
    validateParameters_(parameters_);
  }

  void IDFreeMassErrorEstimator::validateParameters_(const Parameters& p)
  {
    auto positive = [](double value) { return std::isfinite(value) && value > 0.0; };
    if (!positive(p.rt_window_seconds) || !positive(p.precursor_candidate_ppm) || !positive(p.fragment_match_da))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "RT, precursor and fragment matching windows must be finite and positive.");
    }
    if (!positive(p.fragment_low_res_match_da) || p.fragment_low_res_match_da <= p.fragment_match_da ||
        !positive(p.fragment_very_low_res_match_da) || p.fragment_very_low_res_match_da <= p.fragment_low_res_match_da)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Fragment matching windows must increase from narrow to low-resolution to very-low-resolution.");
    }
    if (p.top_peaks < 8 || p.min_matched_peaks < 3 || p.max_candidates_per_bin < 1 ||
        !(p.min_overlap_fraction > 0.0 && p.min_overlap_fraction <= 1.0))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Invalid peak-count, candidate-count or overlap threshold.");
    }
    if (p.min_precursor_clusters < 1 || p.min_precursor_cluster_size < 3 ||
        p.min_tolerance_pairs < p.min_spectrum_pairs || p.min_tolerance_clusters < p.min_precursor_clusters ||
        p.min_fragment_tolerance_pairs < p.min_fragment_pairs || p.min_fragment_tolerance_spectra < 1)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Invalid estimator support thresholds.");
    }
    if (!positive(p.tolerance_sigma_multiplier) || !positive(p.fragment_tolerance_sigma_multiplier) ||
        !positive(p.fragment_high_res_max_sigma_da) || !positive(p.fragment_high_res_max_sigma_ppm) ||
        !(std::isfinite(p.fragment_high_intensity_quantile) && p.fragment_high_intensity_quantile > 0.0 &&
          p.fragment_high_intensity_quantile < 1.0) ||
        !positive(p.fragment_low_res_min_sigma_da) || !positive(p.fragment_low_res_min_sigma_ppm))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Sigma multipliers and fragment-regime thresholds must be finite and positive, and the high-intensity quantile must be in (0,1).");
    }
  }

  void IDFreeMassErrorEstimator::setParameters(const Parameters& parameters)
  {
    validateParameters_(parameters);
    parameters_ = parameters;
    reset();
  }

  const IDFreeMassErrorEstimator::Parameters& IDFreeMassErrorEstimator::getParameters() const
  {
    return parameters_;
  }

  void IDFreeMassErrorEstimator::reset()
  {
    precursor_bins_.clear();
    fragment_bins_.clear();
    precursor_errors_ppm_.clear();
    precursor_errors_da_.clear();
    fragment_errors_da_.clear();
    fragment_errors_ppm_.clear();
    fragment_mixture_errors_ppm_.clear();
    fragment_mixture_mean_mz_.clear();
    fragment_mixture_shared_relative_intensity_.clear();
    fragment_low_res_errors_da_.clear();
    fragment_low_res_errors_ppm_.clear();
    fragment_very_low_res_errors_da_.clear();
    fragment_very_low_res_errors_ppm_.clear();
    diagnostics_ = Diagnostics{};
  }

  std::optional<IDFreeMassErrorEstimator::RobustError> IDFreeMassErrorEstimator::robustError_(const std::vector<double>& values)
  {
    auto finite = sortedFinite(values);
    if (finite.size() < 2) return std::nullopt;

    const double center = Math::median(finite.begin(), finite.end(), true);
    std::vector<double> centered_absolute;
    std::vector<double> absolute_values;
    centered_absolute.reserve(finite.size());
    absolute_values.reserve(finite.size());
    for (double value : finite)
    {
      centered_absolute.push_back(std::abs(value - center));
      absolute_values.push_back(std::abs(value));
    }
    std::sort(centered_absolute.begin(), centered_absolute.end());
    std::sort(absolute_values.begin(), absolute_values.end());

    const double mad = Math::median(centered_absolute.begin(), centered_absolute.end(), true);
    double pair_sigma = 1.4826 * mad;
    if (pair_sigma == 0.0 && finite.size() >= 4)
    {
      pair_sigma = (quantileSorted(finite, 0.75) - quantileSorted(finite, 0.25)) / 1.349;
    }

    RobustError result;
    result.pairwise_median = center;
    result.pairwise_sigma = pair_sigma;
    result.single_measurement_sigma = pair_sigma / std::sqrt(2.0);
    result.pairwise_p95_abs = quantileSorted(absolute_values, 0.95);
    result.robust_inlier_threshold_3sigma = 3.0 * pair_sigma;

    std::vector<double> inlier_centered_absolute;
    inlier_centered_absolute.reserve(centered_absolute.size());
    if (pair_sigma > 0.0)
    {
      for (double value : centered_absolute)
      {
        if (value <= result.robust_inlier_threshold_3sigma) inlier_centered_absolute.push_back(value);
      }
    }
    else
    {
      for (double value : centered_absolute)
      {
        if (value == 0.0) inlier_centered_absolute.push_back(value);
      }
    }

    result.robust_inlier_count = inlier_centered_absolute.size();
    result.robust_outlier_count = finite.size() - result.robust_inlier_count;
    result.robust_inlier_fraction = static_cast<double>(result.robust_inlier_count) / static_cast<double>(finite.size());
    result.robust_inlier_p95_abs_centered = inlier_centered_absolute.empty()
      ? std::numeric_limits<double>::quiet_NaN()
      : quantileSorted(inlier_centered_absolute, 0.95);
    return result;
  }

  std::optional<IDFreeMassErrorEstimator::GaussianUniformFit> IDFreeMassErrorEstimator::gaussianUniformFit_(
    const std::vector<double>& values)
  {
    auto finite_values = sortedFinite(values);
    if (finite_values.size() < 2) return std::nullopt;

    const double uniform_range = finite_values.back() - finite_values.front();
    if (!(std::isfinite(uniform_range) && uniform_range > 0.0)) return std::nullopt;
    const double uniform_density = 1.0 / uniform_range;

    const Size zero_count = static_cast<Size>(std::count(finite_values.begin(), finite_values.end(), 0.0));
    if (2 * zero_count >= finite_values.size()) return std::nullopt;

    const double initial_mean = Math::median(finite_values.begin(), finite_values.end(), true);

    std::vector<double> centered_absolute;
    centered_absolute.reserve(finite_values.size());
    for (double value : finite_values) centered_absolute.push_back(std::abs(value - initial_mean));
    std::sort(centered_absolute.begin(), centered_absolute.end());

    double initial_sigma = 1.4826 * Math::median(centered_absolute.begin(), centered_absolute.end(), true);
    if (!(std::isfinite(initial_sigma) && initial_sigma > 0.0) && finite_values.size() >= 4)
    {
      initial_sigma = (quantileSorted(finite_values, 0.75) - quantileSorted(finite_values, 0.25)) / 1.349;
    }
    if (!(std::isfinite(initial_sigma) && initial_sigma > 0.0)) return std::nullopt;

    constexpr Size max_iterations = 200;
    constexpr double convergence_tolerance = 1e-8;
    constexpr double min_sigma = 1e-6;
    constexpr double min_probability = 1e-6;
    constexpr double inv_sqrt_two_pi = 0.39894228040143267794;

    auto fit_from_sigma = [&](double seed_sigma) -> std::optional<GaussianUniformFit>
    {
      double mean = initial_mean;
      double sigma = std::max(seed_sigma, min_sigma);
      double signal_fraction = 0.5;
      double previous_log_likelihood = -std::numeric_limits<double>::infinity();
      GaussianUniformFit fit;

      for (Size iteration = 1; iteration <= max_iterations; ++iteration)
      {
        double responsibility_sum = 0.0;
        double weighted_value_sum = 0.0;
        double log_likelihood = 0.0;
        std::vector<double> responsibilities(finite_values.size(), 0.0);

        for (Size i = 0; i < finite_values.size(); ++i)
        {
          const double z = (finite_values[i] - mean) / sigma;
          const double gaussian_density = inv_sqrt_two_pi / sigma * std::exp(-0.5 * z * z);
          const double gaussian_term = signal_fraction * gaussian_density;
          const double uniform_term = (1.0 - signal_fraction) * uniform_density;
          const double total_density = gaussian_term + uniform_term;
          if (!(std::isfinite(total_density) && total_density > 0.0)) return std::nullopt;

          const double responsibility = gaussian_term / total_density;
          responsibilities[i] = responsibility;
          responsibility_sum += responsibility;
          weighted_value_sum += responsibility * finite_values[i];
          log_likelihood += std::log(total_density);
        }

        if (!(std::isfinite(responsibility_sum) && responsibility_sum > 0.0)) return std::nullopt;
        const double updated_mean = weighted_value_sum / responsibility_sum;

        double weighted_variance_sum = 0.0;
        for (Size i = 0; i < finite_values.size(); ++i)
        {
          const double centered = finite_values[i] - updated_mean;
          weighted_variance_sum += responsibilities[i] * centered * centered;
        }
        const double updated_sigma = std::sqrt(std::max(weighted_variance_sum / responsibility_sum, min_sigma * min_sigma));
        const double updated_signal_fraction = std::clamp(
          responsibility_sum / static_cast<double>(finite_values.size()), min_probability, 1.0 - min_probability);
        if (!(std::isfinite(updated_mean) && std::isfinite(updated_sigma) && std::isfinite(updated_signal_fraction)))
        {
          return std::nullopt;
        }

        fit.mean = updated_mean;
        fit.pairwise_sigma = updated_sigma;
        fit.signal_fraction = updated_signal_fraction;
        fit.iterations = iteration;
        fit.log_likelihood = log_likelihood;

        if (iteration > 1 &&
            std::abs(log_likelihood - previous_log_likelihood) <=
              convergence_tolerance * (1.0 + std::abs(previous_log_likelihood)))
        {
          fit.converged = true;
          return fit;
        }

        mean = updated_mean;
        sigma = updated_sigma;
        signal_fraction = updated_signal_fraction;
        previous_log_likelihood = log_likelihood;
      }
      return fit;
    };

    std::optional<GaussianUniformFit> best_fit;
    for (double scale : {1.0, 0.5, 0.25, 0.125})
    {
      auto fit = fit_from_sigma(initial_sigma * scale);
      if (fit && (!best_fit || fit->log_likelihood > best_fit->log_likelihood)) best_fit = std::move(fit);
    }
    if (!best_fit || !best_fit->converged || !(best_fit->pairwise_sigma > 0.0)) return std::nullopt;
    return best_fit;
  }

  IDFreeMassErrorEstimator::RobustError IDFreeMassErrorEstimator::errorFromGaussianUniformFit_(
    const std::vector<double>& values,
    const GaussianUniformFit& fit)
  {
    auto finite = sortedFinite(values);
    RobustError result;
    result.pairwise_median = Math::median(finite.begin(), finite.end(), true);
    result.pairwise_sigma = fit.pairwise_sigma;
    result.single_measurement_sigma = fit.pairwise_sigma / std::sqrt(2.0);

    std::vector<double> absolute_values;
    std::vector<double> centered_absolute;
    absolute_values.reserve(finite.size());
    centered_absolute.reserve(finite.size());
    for (double value : finite)
    {
      absolute_values.push_back(std::abs(value));
      centered_absolute.push_back(std::abs(value - fit.mean));
    }
    std::sort(absolute_values.begin(), absolute_values.end());
    std::sort(centered_absolute.begin(), centered_absolute.end());
    result.pairwise_p95_abs = quantileSorted(absolute_values, 0.95);
    result.robust_inlier_threshold_3sigma = 3.0 * fit.pairwise_sigma;

    std::vector<double> inlier_centered_absolute;
    inlier_centered_absolute.reserve(centered_absolute.size());
    for (double value : centered_absolute)
    {
      if (value <= result.robust_inlier_threshold_3sigma) inlier_centered_absolute.push_back(value);
    }
    result.robust_inlier_count = inlier_centered_absolute.size();
    result.robust_outlier_count = finite.size() - result.robust_inlier_count;
    result.robust_inlier_fraction = finite.empty()
      ? 0.0
      : static_cast<double>(result.robust_inlier_count) / static_cast<double>(finite.size());
    result.robust_inlier_p95_abs_centered = inlier_centered_absolute.empty()
      ? std::numeric_limits<double>::quiet_NaN()
      : quantileSorted(inlier_centered_absolute, 0.95);
    return result;
  }

  std::pair<std::vector<double>, std::vector<double>> IDFreeMassErrorEstimator::centroidTopPeaks_(const MSSpectrum& spectrum, Size limit)
  {
    std::vector<PeakValue> peaks;
    peaks.reserve(spectrum.size());
    double base_peak_intensity = 0.0;
    for (const auto& peak : spectrum)
    {
      const double mz = peak.getMZ();
      const double intensity = peak.getIntensity();
      if (std::isfinite(mz) && mz > 0.0 && std::isfinite(intensity) && intensity > 0.0)
      {
        peaks.push_back({mz, intensity});
        base_peak_intensity = std::max(base_peak_intensity, intensity);
      }
    }
    if (peaks.size() > limit)
    {
      std::partial_sort(peaks.begin(), peaks.begin() + limit, peaks.end(), [](const PeakValue& a, const PeakValue& b) {
        if (a.intensity != b.intensity) return a.intensity > b.intensity;
        return a.mz < b.mz;
      });
      peaks.resize(limit);
    }
    std::sort(peaks.begin(), peaks.end(), [](const PeakValue& a, const PeakValue& b) { return a.mz < b.mz; });

    std::vector<double> mz;
    std::vector<double> relative_intensity;
    mz.reserve(peaks.size());
    relative_intensity.reserve(peaks.size());
    for (const auto& peak : peaks)
    {
      mz.push_back(peak.mz);
      relative_intensity.push_back(base_peak_intensity > 0.0 ? peak.intensity / base_peak_intensity : 0.0);
    }
    return {std::move(mz), std::move(relative_intensity)};
  }

  std::pair<std::vector<double>, std::vector<double>> IDFreeMassErrorEstimator::profilePeakCenters_(const MSSpectrum& spectrum, Size limit)
  {
    std::vector<PeakValue> peaks;
    peaks.reserve(spectrum.size());
    double base_peak_intensity = 0.0;
    for (const auto& peak : spectrum)
    {
      const double mz = peak.getMZ();
      const double intensity = peak.getIntensity();
      if (std::isfinite(mz) && mz > 0.0 && std::isfinite(intensity) && intensity >= 0.0)
      {
        peaks.push_back({mz, intensity});
        base_peak_intensity = std::max(base_peak_intensity, intensity);
      }
    }
    if (peaks.size() < 3) return {};
    if (!std::is_sorted(peaks.begin(), peaks.end(), [](const PeakValue& a, const PeakValue& b) { return a.mz < b.mz; }))
    {
      std::stable_sort(peaks.begin(), peaks.end(), [](const PeakValue& a, const PeakValue& b) { return a.mz < b.mz; });
    }

    std::vector<Size> maxima;
    for (Size i = 1; i + 1 < peaks.size(); ++i)
    {
      if (peaks[i].intensity > peaks[i - 1].intensity && peaks[i].intensity >= peaks[i + 1].intensity && peaks[i].intensity > 0.0)
      {
        maxima.push_back(i);
      }
    }
    if (maxima.size() > limit)
    {
      std::partial_sort(maxima.begin(), maxima.begin() + limit, maxima.end(), [&peaks](Size a, Size b) {
        if (peaks[a].intensity != peaks[b].intensity) return peaks[a].intensity > peaks[b].intensity;
        return a < b;
      });
      maxima.resize(limit);
    }

    struct CenterValue
    {
      double mz{0.0};
      double intensity{0.0};
    };
    std::vector<CenterValue> centers;
    centers.reserve(maxima.size());
    constexpr double min_neighbor_fraction = 0.01;
    for (Size i : maxima)
    {
      const double apex = peaks[i].intensity;
      const double left_intensity = peaks[i - 1].intensity;
      const double right_intensity = peaks[i + 1].intensity;
      if (left_intensity <= 0.0 || right_intensity <= 0.0 ||
          left_intensity < apex * min_neighbor_fraction || right_intensity < apex * min_neighbor_fraction)
      {
        continue;
      }

      const double center_sample = peaks[i].mz;
      const double x_left = peaks[i - 1].mz - center_sample;
      const double x_right = peaks[i + 1].mz - center_sample;
      const double y_left = std::log(left_intensity) - std::log(apex);
      const double y_right = std::log(right_intensity) - std::log(apex);
      const double determinant = x_left * x_right * (x_left - x_right);
      if (determinant == 0.0) continue;

      const double a = (y_left * x_right - y_right * x_left) / determinant;
      const double b = (x_left * x_left * y_right - x_right * x_right * y_left) / determinant;
      if (!std::isfinite(a) || !std::isfinite(b) || a >= 0.0) continue;
      const double center_offset = -b / (2.0 * a);
      const double center = center_sample + center_offset;
      const double fitted_log_intensity = std::log(apex) + a * center_offset * center_offset + b * center_offset;
      const double fitted_intensity = std::exp(fitted_log_intensity);
      if (std::isfinite(center) && center >= peaks[i - 1].mz && center <= peaks[i + 1].mz &&
          std::isfinite(fitted_intensity) && fitted_intensity > 0.0)
      {
        centers.push_back({center, fitted_intensity});
      }
    }
    std::sort(centers.begin(), centers.end(), [](const CenterValue& a, const CenterValue& b) { return a.mz < b.mz; });

    std::vector<double> mz;
    std::vector<double> relative_intensity;
    mz.reserve(centers.size());
    relative_intensity.reserve(centers.size());
    for (const auto& center : centers)
    {
      mz.push_back(center.mz);
      relative_intensity.push_back(base_peak_intensity > 0.0 ? center.intensity / base_peak_intensity : 0.0);
    }
    return {std::move(mz), std::move(relative_intensity)};
  }

  std::pair<std::vector<double>, std::vector<double>> IDFreeMassErrorEstimator::matchedDeltas_(
    const std::vector<double>& left_mz,
    const std::vector<double>& right_mz,
    double tolerance_da)
  {
    Size i = 0;
    Size j = 0;
    std::vector<double> deltas;
    std::vector<double> means;
    while (i < left_mz.size() && j < right_mz.size())
    {
      double delta = right_mz[j] - left_mz[i];
      if (std::abs(delta) <= tolerance_da)
      {
        if (j + 1 < right_mz.size())
        {
          const double next_delta = right_mz[j + 1] - left_mz[i];
          if (std::abs(next_delta) < std::abs(delta) && std::abs(next_delta) <= tolerance_da)
          {
            ++j;
            continue;
          }
        }
        if (i + 1 < left_mz.size())
        {
          const double next_delta = right_mz[j] - left_mz[i + 1];
          if (std::abs(next_delta) < std::abs(delta) && std::abs(next_delta) <= tolerance_da)
          {
            ++i;
            continue;
          }
        }
        deltas.push_back(delta);
        means.push_back((left_mz[i] + right_mz[j]) / 2.0);
        ++i;
        ++j;
      }
      else if (delta < -tolerance_da)
      {
        ++j;
      }
      else
      {
        ++i;
      }
    }
    return {std::move(deltas), std::move(means)};
  }

  std::pair<std::vector<double>, std::vector<double>> IDFreeMassErrorEstimator::coarseBinnedDeltas_(
    const std::vector<double>& left_mz,
    const std::vector<double>& left_relative_intensity,
    const std::vector<double>& right_mz,
    const std::vector<double>& right_relative_intensity,
    double bin_width_da,
    std::vector<double>& shared_relative_intensity)
  {

    shared_relative_intensity.clear();
    if (left_relative_intensity.size() != left_mz.size() || right_relative_intensity.size() != right_mz.size()) return {};
    Size i = 0;
    Size j = 0;
    std::vector<double> deltas;
    std::vector<double> means;
    while (i < left_mz.size() && j < right_mz.size())
    {
      const auto left_bin = static_cast<Int64>(std::floor(left_mz[i] / bin_width_da));
      const auto right_bin = static_cast<Int64>(std::floor(right_mz[j] / bin_width_da));
      if (left_bin < right_bin)
      {
        ++i;
        continue;
      }
      if (right_bin < left_bin)
      {
        ++j;
        continue;
      }

      Size next_i = i + 1;
      while (next_i < left_mz.size() &&
             static_cast<Int64>(std::floor(left_mz[next_i] / bin_width_da)) == left_bin)
      {
        ++next_i;
      }
      Size next_j = j + 1;
      while (next_j < right_mz.size() &&
             static_cast<Int64>(std::floor(right_mz[next_j] / bin_width_da)) == right_bin)
      {
        ++next_j;
      }

      if (next_i == i + 1 && next_j == j + 1)
      {
        deltas.push_back(right_mz[j] - left_mz[i]);
        means.push_back((left_mz[i] + right_mz[j]) / 2.0);
        shared_relative_intensity.push_back(std::min(left_relative_intensity[i], right_relative_intensity[j]));
      }
      i = next_i;
      j = next_j;
    }
    return {std::move(deltas), std::move(means)};
  }

  void IDFreeMassErrorEstimator::commitPrecursorCluster_(PrecursorCluster& cluster)
  {
    if (cluster.committed || cluster.count < parameters_.min_precursor_cluster_size) return;
    precursor_errors_da_.insert(precursor_errors_da_.end(), cluster.pending_da.begin(), cluster.pending_da.end());
    precursor_errors_ppm_.insert(precursor_errors_ppm_.end(), cluster.pending_ppm.begin(), cluster.pending_ppm.end());
    diagnostics_.precursor_paired_spectra += cluster.pending_da.size();
    ++diagnostics_.precursor_clusters_used;
    cluster.pending_da.clear();
    cluster.pending_ppm.clear();
    cluster.committed = true;
  }

  void IDFreeMassErrorEstimator::consumePrecursor_(const MSSpectrum& spectrum, double precursor_mz, Int charge)
  {
    ++diagnostics_.precursor_eligible_ms2;
    const Int coarse = static_cast<Int>(std::floor(precursor_mz));
    PrecursorCluster* best = nullptr;
    double best_score = std::numeric_limits<double>::infinity();

    for (Int bin_id = coarse - 1; bin_id <= coarse + 1; ++bin_id)
    {
      auto& clusters = precursor_bins_[{charge, bin_id}];
      std::erase_if(clusters, [&](const PrecursorCluster& cluster) {
        return std::abs(spectrum.getRT() - cluster.last_rt) > parameters_.rt_window_seconds;
      });
      for (auto& cluster : clusters)
      {
        const double denominator = (precursor_mz + cluster.mean_mz) / 2.0;
        const double ppm = (precursor_mz - cluster.mean_mz) / denominator * 1e6;
        if (std::abs(ppm) > parameters_.precursor_candidate_ppm) continue;
        const double score = std::abs(ppm);
        if (score < best_score)
        {
          best_score = score;
          best = &cluster;
        }
      }
    }

    if (best == nullptr)
    {
      auto& clusters = precursor_bins_[{charge, coarse}];
      clusters.push_back({spectrum.getRT(), precursor_mz, 1, precursor_mz, {}, {}, false});
      if (clusters.size() > parameters_.max_candidates_per_bin)
      {
        clusters.erase(clusters.begin(), clusters.begin() + static_cast<std::ptrdiff_t>(clusters.size() - parameters_.max_candidates_per_bin));
      }
      return;
    }

    const double denominator = (precursor_mz + best->last_mz) / 2.0;
    const double delta_da = precursor_mz - best->last_mz;
    const double delta_ppm = delta_da / denominator * 1e6;
    if (best->committed)
    {
      precursor_errors_da_.push_back(delta_da);
      precursor_errors_ppm_.push_back(delta_ppm);
      ++diagnostics_.precursor_paired_spectra;
    }
    else
    {
      best->pending_da.push_back(delta_da);
      best->pending_ppm.push_back(delta_ppm);
    }

    ++best->count;
    best->mean_mz += (precursor_mz - best->mean_mz) / static_cast<double>(best->count);
    best->last_mz = precursor_mz;
    best->last_rt = spectrum.getRT();
    commitPrecursorCluster_(*best);
  }

  std::pair<std::vector<double>, std::vector<double>> IDFreeMassErrorEstimator::fragmentPeaks_(const MSSpectrum& spectrum)
  {
    auto representation = spectrum.getType();
    if (representation == SpectrumSettings::SpectrumType::UNKNOWN && spectrum.size() > 10)
    {
      representation = PeakTypeEstimator::estimateType(spectrum.begin(), spectrum.end());
    }

    if (representation == SpectrumSettings::SpectrumType::CENTROID)
    {
      ++diagnostics_.fragment_centroid_spectra;
      return centroidTopPeaks_(spectrum, parameters_.top_peaks);
    }
    if (representation == SpectrumSettings::SpectrumType::PROFILE)
    {
      ++diagnostics_.fragment_profile_spectra;
      auto peaks = profilePeakCenters_(spectrum, parameters_.top_peaks);
      diagnostics_.fragment_profile_centroids += peaks.first.size();
      if (peaks.first.size() < parameters_.min_matched_peaks) ++diagnostics_.excluded_profile_peak_pick_failure;
      return peaks;
    }

    ++diagnostics_.excluded_fragment_unknown;
    ++diagnostics_.excluded_fragment_profile_or_unknown;
    return {};
  }

  void IDFreeMassErrorEstimator::consumeFragments_(const MSSpectrum& spectrum, double precursor_mz, Int charge)
  {
    auto [mz, relative_intensity] = fragmentPeaks_(spectrum);
    if (mz.size() < parameters_.min_matched_peaks)
    {
      ++diagnostics_.excluded_low_fragment_peaks;
      return;
    }
    ++diagnostics_.fragment_eligible_ms2;

    const Int coarse = static_cast<Int>(std::floor(precursor_mz));
    const FragmentFingerprint* best_candidate = nullptr;
    std::vector<double> best_deltas;
    std::vector<double> best_means;
    double best_score = -1.0;

    for (Int bin_id = coarse - 1; bin_id <= coarse + 1; ++bin_id)
    {
      auto& queue = fragment_bins_[{charge, bin_id}];
      std::erase_if(queue, [&](const FragmentFingerprint& candidate) {
        return std::abs(spectrum.getRT() - candidate.rt) > parameters_.rt_window_seconds;
      });
      for (const auto& candidate : queue)
      {
        const double denominator = (precursor_mz + candidate.precursor_mz) / 2.0;
        const double precursor_ppm = (precursor_mz - candidate.precursor_mz) / denominator * 1e6;
        if (std::abs(precursor_ppm) > parameters_.precursor_candidate_ppm) continue;

        auto [deltas, means] = matchedDeltas_(candidate.mz, mz, parameters_.fragment_match_da);
        const Size required = std::max(parameters_.min_matched_peaks,
          static_cast<Size>(std::ceil(static_cast<double>(std::min(candidate.mz.size(), mz.size())) * parameters_.min_overlap_fraction)));
        if (deltas.size() < required) continue;

        const double score = static_cast<double>(deltas.size()) / static_cast<double>(std::max(candidate.mz.size(), mz.size()));
        if (score > best_score)
        {
          best_score = score;
          best_candidate = &candidate;
          best_deltas = std::move(deltas);
          best_means = std::move(means);
        }
      }
    }

    if (best_candidate != nullptr)
    {
      fragment_errors_da_.insert(fragment_errors_da_.end(), best_deltas.begin(), best_deltas.end());
      for (Size i = 0; i < best_deltas.size(); ++i)
      {
        if (best_means[i] > 0.0) fragment_errors_ppm_.push_back(best_deltas[i] / best_means[i] * 1e6);
      }

      constexpr double mixture_bin_width_da = 1.0005079;
      std::vector<double> mixture_shared_relative_intensity;
      auto [mixture_deltas, mixture_means] = coarseBinnedDeltas_(
        best_candidate->mz,
        best_candidate->relative_intensity,
        mz,
        relative_intensity,
        mixture_bin_width_da,
        mixture_shared_relative_intensity);
      for (Size i = 0; i < mixture_deltas.size(); ++i)
      {
        if (mixture_means[i] > 0.0)
        {
          fragment_mixture_errors_ppm_.push_back(mixture_deltas[i] / mixture_means[i] * 1e6);
          fragment_mixture_mean_mz_.push_back(mixture_means[i]);
          fragment_mixture_shared_relative_intensity_.push_back(mixture_shared_relative_intensity[i]);
        }
      }

      auto [wide_deltas, wide_means] = matchedDeltas_(best_candidate->mz, mz, parameters_.fragment_low_res_match_da);
      fragment_low_res_errors_da_.insert(fragment_low_res_errors_da_.end(), wide_deltas.begin(), wide_deltas.end());
      for (Size i = 0; i < wide_deltas.size(); ++i)
      {
        if (wide_means[i] > 0.0) fragment_low_res_errors_ppm_.push_back(wide_deltas[i] / wide_means[i] * 1e6);
      }

      auto [very_wide_deltas, very_wide_means] = matchedDeltas_(best_candidate->mz, mz, parameters_.fragment_very_low_res_match_da);
      fragment_very_low_res_errors_da_.insert(fragment_very_low_res_errors_da_.end(), very_wide_deltas.begin(), very_wide_deltas.end());
      for (Size i = 0; i < very_wide_deltas.size(); ++i)
      {
        if (very_wide_means[i] > 0.0) fragment_very_low_res_errors_ppm_.push_back(very_wide_deltas[i] / very_wide_means[i] * 1e6);
      }
      ++diagnostics_.fragment_paired_spectra;
    }

    auto& queue = fragment_bins_[{charge, coarse}];
    queue.push_back({spectrum.getRT(), precursor_mz, std::move(mz), std::move(relative_intensity)});
    while (queue.size() > parameters_.max_candidates_per_bin) queue.pop_front();
  }

  void IDFreeMassErrorEstimator::consumeSpectrum(const MSSpectrum& spectrum)
  {
    if (spectrum.getMSLevel() != 2) return;
    const auto& precursors = spectrum.getPrecursors();
    if (precursors.empty())
    {
      ++diagnostics_.excluded_missing_precursor;
      return;
    }

    const double precursor_mz = precursors.front().getMZ();
    const Int charge = precursors.front().getCharge();
    if (!std::isfinite(precursor_mz) || precursor_mz <= 0.0)
    {
      ++diagnostics_.excluded_missing_precursor;
      return;
    }

    // Precursor precision is charge-state specific. Fragment precision only needs
    // a repeatable precursor/isolation m/z anchor, so charge-less DIA spectra can
    // still contribute fragment evidence without inventing a precursor charge.
    if (charge > 0)
    {
      consumePrecursor_(spectrum, precursor_mz, charge);
    }
    consumeFragments_(spectrum, precursor_mz, charge);
  }

  IDFreeMassErrorEstimator::Result IDFreeMassErrorEstimator::compute(const MSExperiment& experiment)
  {
    reset();
    for (const auto& spectrum : experiment) consumeSpectrum(spectrum);
    return getResult();
  }

  std::optional<IDFreeMassErrorEstimator::RobustError> IDFreeMassErrorEstimator::getPrecursorPrecisionPPM() const
  {
    return robustError_(precursor_errors_ppm_);
  }

  IDFreeMassErrorEstimator::Result IDFreeMassErrorEstimator::getResult() const
  {
    Result result;
    auto precursor_ppm = robustError_(precursor_errors_ppm_);
    auto precursor_da = robustError_(precursor_errors_da_);
    auto fragment_da_narrow = robustError_(fragment_errors_da_);
    auto fragment_ppm_narrow = robustError_(fragment_errors_ppm_);
    auto fragment_ppm_mixture_fit = gaussianUniformFit_(fragment_mixture_errors_ppm_);
    std::optional<RobustError> fragment_ppm_mixture;
    if (fragment_ppm_mixture_fit)
    {
      fragment_ppm_mixture = errorFromGaussianUniformFit_(fragment_mixture_errors_ppm_, *fragment_ppm_mixture_fit);
    }

    std::vector<double> fragment_high_intensity_errors_ppm;
    std::vector<double> fragment_high_intensity_mean_mz;
    double fragment_high_intensity_threshold = 0.0;
    if (!fragment_mixture_shared_relative_intensity_.empty() &&
        fragment_mixture_shared_relative_intensity_.size() == fragment_mixture_errors_ppm_.size() &&
        fragment_mixture_mean_mz_.size() == fragment_mixture_errors_ppm_.size())
    {
      auto finite_relative_intensity = sortedFinite(fragment_mixture_shared_relative_intensity_);
      if (!finite_relative_intensity.empty())
      {
        fragment_high_intensity_threshold = quantileSorted(finite_relative_intensity, parameters_.fragment_high_intensity_quantile);
        for (Size i = 0; i < fragment_mixture_errors_ppm_.size(); ++i)
        {
          if (std::isfinite(fragment_mixture_shared_relative_intensity_[i]) &&
              fragment_mixture_shared_relative_intensity_[i] > fragment_high_intensity_threshold)
          {
            fragment_high_intensity_errors_ppm.push_back(fragment_mixture_errors_ppm_[i]);
            fragment_high_intensity_mean_mz.push_back(fragment_mixture_mean_mz_[i]);
          }
        }
      }
    }
    auto fragment_high_intensity_fit = gaussianUniformFit_(fragment_high_intensity_errors_ppm);
    std::optional<RobustError> fragment_high_intensity_ppm;
    if (fragment_high_intensity_fit)
    {
      fragment_high_intensity_ppm = errorFromGaussianUniformFit_(fragment_high_intensity_errors_ppm, *fragment_high_intensity_fit);
    }

    auto fragment_da_wide = robustError_(fragment_low_res_errors_da_);
    auto fragment_ppm_wide = robustError_(fragment_low_res_errors_ppm_);
    auto fragment_da_very_wide = robustError_(fragment_very_low_res_errors_da_);
    auto fragment_ppm_very_wide = robustError_(fragment_very_low_res_errors_ppm_);

    result.fragment_match_window_da = parameters_.fragment_match_da;
    result.diagnostics = diagnostics_;
    result.diagnostics.fragment_pairs = fragment_errors_da_.size();
    result.diagnostics.fragment_mixture_pairs = fragment_mixture_errors_ppm_.size();
    result.diagnostics.fragment_low_res_pairs = fragment_low_res_errors_da_.size();
    result.diagnostics.fragment_very_low_res_pairs = fragment_very_low_res_errors_da_.size();
    result.diagnostics.fragment_zero_delta_pairs = static_cast<Size>(
      std::count(fragment_mixture_errors_ppm_.begin(), fragment_mixture_errors_ppm_.end(), 0.0));
    if (!fragment_mixture_errors_ppm_.empty())
    {
      result.diagnostics.fragment_zero_delta_fraction =
        static_cast<double>(result.diagnostics.fragment_zero_delta_pairs) / static_cast<double>(fragment_mixture_errors_ppm_.size());
    }
    result.diagnostics.fragment_mixture_rejected_zero_quantization =
      !fragment_mixture_errors_ppm_.empty() &&
      2 * result.diagnostics.fragment_zero_delta_pairs >= fragment_mixture_errors_ppm_.size();
    if (!fragment_mixture_mean_mz_.empty())
    {
      auto reference_mz = sortedFinite(fragment_mixture_mean_mz_);
      result.diagnostics.fragment_mixture_reference_mz = Math::median(reference_mz.begin(), reference_mz.end(), true);
    }
    if (fragment_ppm_mixture_fit)
    {
      result.diagnostics.fragment_mixture_mean_ppm = fragment_ppm_mixture_fit->mean;
      result.diagnostics.fragment_mixture_signal_fraction = fragment_ppm_mixture_fit->signal_fraction;
      result.diagnostics.fragment_mixture_signal_pairs = static_cast<Size>(std::llround(
        fragment_ppm_mixture_fit->signal_fraction * static_cast<double>(fragment_mixture_errors_ppm_.size())));
      result.diagnostics.fragment_mixture_iterations = fragment_ppm_mixture_fit->iterations;
      result.diagnostics.fragment_mixture_converged = fragment_ppm_mixture_fit->converged;
    }

    result.diagnostics.fragment_high_intensity_pairs = fragment_high_intensity_errors_ppm.size();
    result.diagnostics.fragment_high_intensity_threshold = fragment_high_intensity_threshold;
    if (fragment_high_intensity_fit)
    {
      result.diagnostics.fragment_high_intensity_signal_fraction = fragment_high_intensity_fit->signal_fraction;
      result.diagnostics.fragment_high_intensity_signal_pairs = static_cast<Size>(std::llround(
        fragment_high_intensity_fit->signal_fraction * static_cast<double>(fragment_high_intensity_errors_ppm.size())));
      result.diagnostics.fragment_high_intensity_mixture_converged = fragment_high_intensity_fit->converged;
    }

    const RobustError* selected_fragment_da = fragment_da_narrow ? &*fragment_da_narrow : nullptr;
    const RobustError* selected_fragment_ppm = fragment_ppm_narrow ? &*fragment_ppm_narrow : nullptr;
    Size selected_fragment_pairs = fragment_errors_da_.size();

    const double fragment_mixture_single_sigma_da = fragment_ppm_mixture
      ? fragment_ppm_mixture->single_measurement_sigma * result.diagnostics.fragment_mixture_reference_mz * 1e-6
      : 0.0;
    const bool full_mixture_high_resolution = fragment_ppm_mixture &&
      fragment_ppm_mixture->single_measurement_sigma > 0.0 &&
      fragment_ppm_mixture->single_measurement_sigma <= parameters_.fragment_high_res_max_sigma_ppm &&
      fragment_mixture_single_sigma_da <= parameters_.fragment_high_res_max_sigma_da;

    double fragment_high_intensity_reference_mz = 0.0;
    if (!fragment_high_intensity_mean_mz.empty())
    {
      auto reference_mz = sortedFinite(fragment_high_intensity_mean_mz);
      fragment_high_intensity_reference_mz = Math::median(reference_mz.begin(), reference_mz.end(), true);
    }
    const double fragment_high_intensity_single_sigma_da = fragment_high_intensity_ppm
      ? fragment_high_intensity_ppm->single_measurement_sigma * fragment_high_intensity_reference_mz * 1e-6
      : 0.0;
    const bool high_intensity_high_resolution = fragment_high_intensity_ppm &&
      fragment_high_intensity_ppm->single_measurement_sigma > 0.0 &&
      fragment_high_intensity_ppm->single_measurement_sigma <= parameters_.fragment_high_res_max_sigma_ppm &&
      fragment_high_intensity_single_sigma_da <= parameters_.fragment_high_res_max_sigma_da &&
      result.diagnostics.fragment_high_intensity_signal_pairs >= parameters_.min_fragment_pairs;

    if (full_mixture_high_resolution)
    {
      result.fragment_resolution_regime = FragmentResolutionRegime::HIGH_RESOLUTION;
      selected_fragment_ppm = &*fragment_ppm_mixture;
      selected_fragment_pairs = result.diagnostics.fragment_mixture_signal_pairs;
    }
    else if (!result.diagnostics.fragment_mixture_rejected_zero_quantization && high_intensity_high_resolution)
    {
      result.fragment_resolution_regime = FragmentResolutionRegime::HIGH_RESOLUTION;
      selected_fragment_ppm = &*fragment_high_intensity_ppm;
      selected_fragment_pairs = result.diagnostics.fragment_high_intensity_signal_pairs;
      result.diagnostics.fragment_high_intensity_fallback_used = true;
    }
    else if (fragment_da_narrow && fragment_ppm_narrow &&
             fragment_da_narrow->single_measurement_sigma >= parameters_.fragment_low_res_min_sigma_da &&
             fragment_ppm_narrow->single_measurement_sigma >= parameters_.fragment_low_res_min_sigma_ppm)
    {
      result.fragment_resolution_regime = FragmentResolutionRegime::LOW_RESOLUTION;
      struct Candidate
      {
        double window;
        const RobustError* da;
        const RobustError* ppm;
        Size pair_count;
      };
      const Candidate candidates[] = {
          {parameters_.fragment_low_res_match_da,
           fragment_da_wide ? &*fragment_da_wide : nullptr,
           fragment_ppm_wide ? &*fragment_ppm_wide : nullptr,
           fragment_low_res_errors_da_.size()},
          {parameters_.fragment_very_low_res_match_da,
           fragment_da_very_wide ? &*fragment_da_very_wide : nullptr,
           fragment_ppm_very_wide ? &*fragment_ppm_very_wide : nullptr,
           fragment_very_low_res_errors_da_.size()}
        };
      const Candidate* selected = nullptr;
      const Candidate* fallback = nullptr;
      for (const auto& candidate : candidates)
      {
        if (candidate.da == nullptr || candidate.ppm == nullptr) continue;
        fallback = &candidate;
        if (candidate.da->robust_inlier_threshold_3sigma < 0.9 * candidate.window)
        {
          selected = &candidate;
          break;
        }
      }
      if (selected == nullptr) selected = fallback;
      if (selected != nullptr)
      {
        result.fragment_match_window_da = selected->window;
        selected_fragment_da = selected->da;
        selected_fragment_ppm = selected->ppm;
        selected_fragment_pairs = selected->pair_count;
        result.fragment_window_censored = selected->da->robust_inlier_threshold_3sigma >= 0.9 * selected->window;
      }
      else
      {
        result.fragment_window_censored = true;
      }
    }

    const bool enough_precursor = diagnostics_.precursor_paired_spectra >= parameters_.min_spectrum_pairs &&
                                  diagnostics_.precursor_clusters_used >= parameters_.min_precursor_clusters;
    if (enough_precursor && precursor_ppm && precursor_ppm->single_measurement_sigma > 0.0) result.precursor_ppm = precursor_ppm;
    if (enough_precursor && precursor_da && precursor_da->single_measurement_sigma > 0.0) result.precursor_da = precursor_da;

    const bool enough_fragment = fragment_errors_da_.size() >= parameters_.min_fragment_pairs;
    const bool enough_selected_fragment = result.fragment_resolution_regime == FragmentResolutionRegime::HIGH_RESOLUTION
      ? selected_fragment_pairs >= parameters_.min_fragment_pairs
      : enough_fragment;
    if (enough_fragment && selected_fragment_da != nullptr && selected_fragment_da->single_measurement_sigma > 0.0)
    {
      result.fragment_da = *selected_fragment_da;
    }
    if (enough_selected_fragment && selected_fragment_ppm != nullptr && selected_fragment_ppm->single_measurement_sigma > 0.0)
    {
      result.fragment_ppm = *selected_fragment_ppm;
    }

    const bool enough_tolerance_support = enough_precursor &&
                                          diagnostics_.precursor_paired_spectra >= parameters_.min_tolerance_pairs &&
                                          diagnostics_.precursor_clusters_used >= parameters_.min_tolerance_clusters;
    if (enough_tolerance_support && precursor_ppm && precursor_ppm->single_measurement_sigma > 0.0)
    {
      ToleranceSuggestion suggestion;
      suggestion.single_measurement_sigma = precursor_ppm->single_measurement_sigma;
      suggestion.sigma_multiplier = parameters_.tolerance_sigma_multiplier;
      suggestion.tolerance = suggestion.sigma_multiplier * suggestion.single_measurement_sigma;
      suggestion.unit = "ppm";
      suggestion.confidence = diagnostics_.precursor_paired_spectra >= 1000 && diagnostics_.precursor_clusters_used >= 250 ? "high" : "moderate";
      suggestion.support = diagnostics_.precursor_paired_spectra;
      result.precursor_tolerance_ppm = std::move(suggestion);
    }

    const bool enough_fragment_tolerance_support = enough_selected_fragment &&
      selected_fragment_pairs >= parameters_.min_fragment_tolerance_pairs &&
      diagnostics_.fragment_paired_spectra >= parameters_.min_fragment_tolerance_spectra;
    if (enough_fragment_tolerance_support && selected_fragment_da != nullptr && selected_fragment_ppm != nullptr &&
        selected_fragment_da->single_measurement_sigma > 0.0 && selected_fragment_ppm->single_measurement_sigma > 0.0)
    {
      ToleranceSuggestion suggestion;
      suggestion.sigma_multiplier = parameters_.fragment_tolerance_sigma_multiplier;
      suggestion.confidence = selected_fragment_pairs >= 10000 && diagnostics_.fragment_paired_spectra >= 500 ? "high" : "moderate";
      suggestion.support = selected_fragment_pairs;
      if (result.fragment_resolution_regime == FragmentResolutionRegime::HIGH_RESOLUTION)
      {
        suggestion.single_measurement_sigma = selected_fragment_ppm->single_measurement_sigma;
        suggestion.tolerance = suggestion.sigma_multiplier * suggestion.single_measurement_sigma;
        suggestion.unit = "ppm";
        result.fragment_tolerance_ppm = std::move(suggestion);
      }
      else if (result.fragment_resolution_regime == FragmentResolutionRegime::LOW_RESOLUTION && !result.fragment_window_censored)
      {
        suggestion.single_measurement_sigma = selected_fragment_da->single_measurement_sigma;
        suggestion.tolerance = suggestion.sigma_multiplier * suggestion.single_measurement_sigma;
        suggestion.unit = "Da";
        result.fragment_tolerance_da = std::move(suggestion);
      }
    }

    return result;
  }
} // namespace OpenMS
