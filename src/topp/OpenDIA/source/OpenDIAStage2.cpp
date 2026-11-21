// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIAStage2.h"
#include "OpenDIAPeakPicker.h"

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace OpenMS::OpenDIA
{
  namespace
  {
    constexpr std::size_t FRAGMENTS = TraceCollector::FRAGMENTS;
    constexpr std::size_t RT_BINS = TraceCollector::RT_BINS;
    constexpr std::size_t IM_BINS = TraceCollector::IM_BINS;

    double pearson(const std::vector<double>& x, const std::vector<double>& y)
    {
      const std::size_t n = x.size();
      if (n < 3) return 0.0;
      const double mx = std::accumulate(x.begin(), x.end(), 0.0) / n;
      const double my = std::accumulate(y.begin(), y.end(), 0.0) / n;
      double sxy = 0.0, sxx = 0.0, syy = 0.0;
      for (std::size_t i = 0; i < n; ++i)
      {
        sxy += (x[i] - mx) * (y[i] - my);
        sxx += (x[i] - mx) * (x[i] - mx);
        syy += (y[i] - my) * (y[i] - my);
      }
      return sxx > 0.0 && syy > 0.0 ? sxy / std::sqrt(sxx * syy) : 0.0;
    }

    /// IM-summed fragment traces of one window over its scan bins, and the peak group picked around the centre bin.
    struct PickedWindow
    {
      std::vector<std::size_t> bins;      ///< RT bins that received an isolating MS2 scan.
      std::vector<std::size_t> channels;  ///< Non-padding fragment channels.
      std::vector<float> traces;          ///< channels x n, fragment-major.
      std::size_t n = 0;
      PickedPeak peak;                    ///< peak.found is false with fewer than 3 scan bins, no channel or no signal.
    };

    PickedWindow pickWindow(const float* tensor, const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence)
    {
      PickedWindow w;
      // Scan bins only: a bin without an isolating scan carries no information, not zero intensity.
      for (std::size_t b = 0; b < RT_BINS; ++b)
      {
        if (evidence.ms2_mask >> b & 1u) w.bins.push_back(b);
      }
      for (std::size_t c = 0; c < FRAGMENTS; ++c)
      {
        if (fragment_mz[c] != 0.0) w.channels.push_back(c);
      }
      w.n = w.bins.size();
      if (w.n < 3 || w.channels.empty()) return w;
      w.traces.assign(w.channels.size() * w.n, 0.0f);
      for (std::size_t k = 0; k < w.channels.size(); ++k)
      {
        for (std::size_t s = 0; s < w.n; ++s)
        {
          float v = 0.0f;
          for (std::size_t m = 0; m < IM_BINS; ++m) v += tensor[(w.channels[k] * RT_BINS + w.bins[s]) * IM_BINS + m];
          w.traces[k * w.n + s] = v;
        }
      }
      // The window is centred on the candidate apex: hint = the scan bin nearest the centre bin.
      std::size_t hint = 0;
      const auto distance = [](std::size_t b) { return b > RT_BINS / 2 ? b - RT_BINS / 2 : RT_BINS / 2 - b; };
      for (std::size_t s = 1; s < w.n; ++s)
      {
        if (distance(w.bins[s]) < distance(w.bins[hint])) hint = s;
      }
      std::vector<double> scratch;
      w.peak = pickPeakGroup(w.traces.data(), w.channels.size(), w.n, w.n, 1, hint, PeakPickerOptions{}, scratch);
      return w;
    }

    /// Solve A x = b for a symmetric positive definite A (Cholesky). A is n x n row-major.
    std::vector<double> solveSPD(std::vector<double> a, std::vector<double> b, std::size_t n)
    {
      for (std::size_t j = 0; j < n; ++j)
      {
        double d = a[j * n + j];
        for (std::size_t k = 0; k < j; ++k) d -= a[j * n + k] * a[j * n + k];
        if (!(d > 0.0)) throw std::invalid_argument("Second-stage covariance is not positive definite.");
        a[j * n + j] = std::sqrt(d);
        for (std::size_t i = j + 1; i < n; ++i)
        {
          double s = a[i * n + j];
          for (std::size_t k = 0; k < j; ++k) s -= a[i * n + k] * a[j * n + k];
          a[i * n + j] = s / a[j * n + j];
        }
      }
      for (std::size_t i = 0; i < n; ++i)
      {
        for (std::size_t k = 0; k < i; ++k) b[i] -= a[i * n + k] * b[k];
        b[i] /= a[i * n + i];
      }
      for (std::size_t i = n; i-- > 0;)
      {
        for (std::size_t k = i + 1; k < n; ++k) b[i] -= a[k * n + i] * b[k];
        b[i] /= a[i * n + i];
      }
      return b;
    }

    /// Indices of precursor winners (highest score, ties to the lower index) among the given rows.
    std::vector<std::size_t> winners(const std::vector<Stage2Row>& rows, const std::vector<std::size_t>& subset,
                                     const std::vector<double>& score)
    {
      std::unordered_map<std::string, std::size_t> best;
      best.reserve(subset.size());
      for (std::size_t i : subset)
      {
        auto [it, inserted] = best.try_emplace(rows[i].precursor_id, i);
        if (!inserted && score[i] > score[it->second]) it->second = i;
      }
      std::vector<std::size_t> out;
      out.reserve(best.size());
      for (const auto& entry : best) out.push_back(entry.second);
      std::sort(out.begin(), out.end());
      return out;
    }

    /// Target winners at decoy-count q = (D+1)/T <= fdr (exact ties grouped), and all decoy winners.
    void selectTraining(const std::vector<Stage2Row>& rows, const std::vector<std::size_t>& winner,
                        const std::vector<double>& score, double fdr,
                        std::vector<std::size_t>& positives, std::vector<std::size_t>& negatives)
    {
      std::vector<std::size_t> order = winner;
      std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b)
      { return score[a] != score[b] ? score[a] > score[b] : a < b; });
      // q at each position, then reverse cumulative minimum.
      std::vector<double> q(order.size(), 1.0);
      std::size_t targets = 0, decoys = 0;
      for (std::size_t k = 0; k < order.size();)
      {
        std::size_t end = k;
        while (end < order.size() && score[order[end]] == score[order[k]])
        {
          (rows[order[end]].is_decoy ? decoys : targets) += 1;
          ++end;
        }
        const double value = targets == 0 ? 1.0 : std::min(1.0, (decoys + 1.0) / targets);
        for (std::size_t m = k; m < end; ++m) q[m] = value;
        k = end;
      }
      for (std::size_t k = order.size(); k-- > 1;) q[k - 1] = std::min(q[k - 1], q[k]);
      positives.clear();
      negatives.clear();
      for (std::size_t k = 0; k < order.size(); ++k)
      {
        if (rows[order[k]].is_decoy) negatives.push_back(order[k]);
        else if (q[k] <= fdr) positives.push_back(order[k]);
      }
    }
  }

  const std::array<std::string, STAGE2_FEATURES>& stage2FeatureNames()
  {
    static const std::array<std::string, STAGE2_FEATURES> names = {
      "cnn_z", "fragment_fraction", "fragment_coelution", "library_correlation",
      "ms2_abs_ppm", "ms1_abs_ppm", "ms1_coelution", "ms1_isotope"};
    return names;
  }

  std::array<double, STAGE2_FEATURES> stage2Features(const float* tensor, const std::array<double, 6>& library_intensities,
                                                     const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence,
                                                     const Stage2ExtractionOptions& options)
  {
    std::array<double, STAGE2_FEATURES> f{};
    const double width = 2.0 * options.rt_half_width / RT_BINS;
    // Core bins that received an MS2 scan isolating the precursor (zeros elsewhere mean no scan, not no signal).
    std::vector<std::size_t> ms2_bins, ms1_bins;
    for (std::size_t b = 0; b < RT_BINS; ++b)
    {
      const double centre = -options.rt_half_width + (b + 0.5) * width;
      if (std::abs(centre) > options.core_half_width) continue;
      if (evidence.ms2_mask >> b & 1u) ms2_bins.push_back(b);
      if (evidence.ms1_mask >> b & 1u) ms1_bins.push_back(b);
    }
    // Fragment traces over the core MS2 bins, summed over ion mobility.
    std::vector<std::vector<double>> traces;
    std::vector<double> areas, library;
    std::vector<double> summed(RT_BINS, 0.0);
    std::size_t valid = 0, with_signal = 0;
    for (std::size_t c = 0; c < FRAGMENTS; ++c)
    {
      if (fragment_mz[c] == 0.0) continue;
      ++valid;
      std::vector<double> trace;
      double area = 0.0;
      for (std::size_t b : ms2_bins)
      {
        double v = 0.0;
        for (std::size_t m = 0; m < IM_BINS; ++m) v += tensor[(c * RT_BINS + b) * IM_BINS + m];
        trace.push_back(v);
        summed[b] += v;
        area += v;
      }
      areas.push_back(area);
      library.push_back(library_intensities[c]);
      if (area > 0.0)
      {
        ++with_signal;
        traces.push_back(std::move(trace));
      }
    }
    f[1] = valid ? static_cast<double>(with_signal) / valid : 0.0;
    if (traces.size() >= 2 && ms2_bins.size() >= 3)
    {
      double sum = 0.0;
      std::size_t pairs = 0;
      for (std::size_t i = 0; i < traces.size(); ++i)
      {
        for (std::size_t j = i + 1; j < traces.size(); ++j)
        {
          sum += pearson(traces[i], traces[j]);
          ++pairs;
        }
      }
      f[2] = sum / pairs;
    }
    if (valid >= 3 && with_signal >= 2)
    {
      std::vector<double> a, l;
      for (std::size_t i = 0; i < areas.size(); ++i)
      {
        a.push_back(std::sqrt(areas[i]));
        l.push_back(std::sqrt(std::max(0.0, library[i])));
      }
      f[3] = pearson(a, l);
    }
    f[4] = evidence.ms2_weight > 0.0
             ? std::min(options.ppm_cap, std::abs(evidence.ms2_ppm_sum / evidence.ms2_weight - options.ms2_offset_ppm))
             : options.ppm_cap;
    f[5] = evidence.ms1_weight > 0.0
             ? std::min(options.ppm_cap, std::abs(evidence.ms1_ppm_sum / evidence.ms1_weight - options.ms1_offset_ppm))
             : options.ppm_cap;
    // MS1 trace against the summed fragment trace, linearly interpolated between the MS2 bins around each MS1 bin.
    if (ms1_bins.size() >= 3 && ms2_bins.size() >= 2)
    {
      std::vector<double> x, y;
      for (std::size_t b : ms1_bins)
      {
        const auto hi = std::lower_bound(ms2_bins.begin(), ms2_bins.end(), b);
        double value;
        if (hi == ms2_bins.begin()) value = summed[*hi];
        else if (hi == ms2_bins.end()) value = summed[ms2_bins.back()];
        else if (*hi == b) value = summed[b];
        else
        {
          const std::size_t lo = *(hi - 1);
          const double t = static_cast<double>(b - lo) / static_cast<double>(*hi - lo);
          value = (1.0 - t) * summed[lo] + t * summed[*hi];
        }
        x.push_back(evidence.ms1[b]);
        y.push_back(value);
      }
      f[6] = pearson(x, y);
    }
    f[7] = evidence.ms1_weight > 0.0 && evidence.ms1_isotope > 0.0f ? 1.0 : 0.0;
    return f;
  }

  Stage2FeatureSet parseStage2FeatureSet(const std::string& name)
  {
    Stage2FeatureSet set{false, false, false, false};
    const std::string spec = name == "extraction_search" ? std::string("extraction+search") : name;
    std::size_t start = 0;
    while (true)
    {
      const std::size_t plus = spec.find('+', start);
      const std::string group = spec.substr(start, plus == std::string::npos ? std::string::npos : plus - start);
      bool* flag = group == "extraction" ? &set.extraction : group == "search" ? &set.search
                   : group == "boundary" ? &set.boundary : group == "byion" ? &set.byion : nullptr;
      if (flag == nullptr || *flag) throw std::invalid_argument("Unknown or repeated second-stage feature group '" + group + "' in " + name);
      *flag = true;
      if (plus == std::string::npos) break;
      start = plus + 1;
    }
    return set;
  }

  const std::array<std::string, STAGE2_SEARCH_FEATURES>& stage2SearchFeatureNames()
  {
    static const std::array<std::string, STAGE2_SEARCH_FEATURES> names = {
      "transition_support", "log_matched_scans", "log_peak_width", "abs_library_rt_error"};
    return names;
  }

  std::vector<std::string> stage2FeatureSetNames(const Stage2FeatureSet& set)
  {
    std::vector<std::string> names;
    if (set.extraction)
    {
      names.assign(stage2FeatureNames().begin(), stage2FeatureNames().end());
    }
    else
    {
      names.push_back(stage2FeatureNames()[0]);
    }
    if (set.search)
    {
      names.insert(names.end(), stage2SearchFeatureNames().begin(), stage2SearchFeatureNames().end());
    }
    if (set.boundary)
    {
      names.insert(names.end(), stage2BoundaryFeatureNames().begin(), stage2BoundaryFeatureNames().end());
    }
    if (set.byion)
    {
      names.insert(names.end(), stage2ByIonFeatureNames().begin(), stage2ByIonFeatureNames().end());
    }
    return names;
  }

  const std::array<std::string, STAGE2_BYION_FEATURES>& stage2ByIonFeatureNames()
  {
    static const std::array<std::string, STAGE2_BYION_FEATURES> names = {"byion_b", "byion_y"};
    return names;
  }

  std::array<double, STAGE2_BYION_FEATURES> stage2ByIonFeatures(const TraceEvidence& evidence)
  {
    return {static_cast<double>(evidence.byion_b), static_cast<double>(evidence.byion_y)};
  }

  const std::array<std::string, STAGE2_BOUNDARY_FEATURES>& stage2BoundaryFeatureNames()
  {
    static const std::array<std::string, STAGE2_BOUNDARY_FEATURES> names = {
      "boundary_log_fwhm", "boundary_log_width", "boundary_coelution", "boundary_library_correlation", "boundary_apex_agreement"};
    return names;
  }

  std::array<double, STAGE2_BOUNDARY_FEATURES> stage2BoundaryFeatures(const float* tensor, const std::array<double, 6>& library_intensities,
                                                                      const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence,
                                                                      const Stage2ExtractionOptions& options)
  {
    std::array<double, STAGE2_BOUNDARY_FEATURES> f{};
    const PickedWindow picked = pickWindow(tensor, fragment_mz, evidence);
    if (!picked.peak.found) return f;
    const std::size_t n = picked.n;
    const auto& channels = picked.channels;
    const auto& traces = picked.traces;
    const auto& peak = picked.peak;
    const double width = 2.0 * options.rt_half_width / RT_BINS;
    // Positions in seconds: bin centres of the scan bins, linearly interpolated at fractional indices.
    std::vector<double> centres(n);
    for (std::size_t s = 0; s < n; ++s) centres[s] = -options.rt_half_width + (picked.bins[s] + 0.5) * width;
    f[0] = std::log1p(std::max(0.0, interpolateAt(centres, peak.half_right) - interpolateAt(centres, peak.half_left)));
    f[1] = std::log1p(std::max(0.0, centres[peak.right] - centres[peak.left]));
    const std::size_t span = peak.right - peak.left + 1;
    std::vector<std::vector<double>> inside;
    std::vector<double> areas, library;
    for (std::size_t k = 0; k < channels.size(); ++k)
    {
      std::vector<double> trace(traces.begin() + k * n + peak.left, traces.begin() + k * n + peak.right + 1);
      const double area = std::accumulate(trace.begin(), trace.end(), 0.0);
      areas.push_back(std::sqrt(area));
      library.push_back(std::sqrt(std::max(0.0, library_intensities[channels[k]])));
      if (area > 0.0) inside.push_back(std::move(trace));
    }
    if (inside.size() >= 2 && span >= 3)
    {
      double sum = 0.0;
      std::size_t pairs = 0;
      for (std::size_t i = 0; i < inside.size(); ++i)
      {
        for (std::size_t j = i + 1; j < inside.size(); ++j)
        {
          sum += pearson(inside[i], inside[j]);
          ++pairs;
        }
      }
      f[2] = sum / pairs;
    }
    if (channels.size() >= 3 && inside.size() >= 2) f[3] = pearson(areas, library);
    f[4] = static_cast<double>(peak.agreeing_fragments) / static_cast<double>(channels.size());
    return f;
  }

  const std::array<std::string, STAGE2_PEAKGROUP_FEATURES>& stage2PeakGroupFeatureNames()
  {
    static const std::array<std::string, STAGE2_PEAKGROUP_FEATURES> names = {
      "pg_xcorr_coelution", "pg_xcorr_coelution_w", "pg_xcorr_shape", "pg_xcorr_shape_w", "pg_log_sn",
      "pg_intensity_fraction", "pg_massdev", "pg_massdev_w", "pg_isotope_corr", "pg_isotope_overlap"};
    return names;
  }

  const std::array<std::string, STAGE2_EXTRA_FEATURES>& stage2ExtraFeatureNames()
  {
    static const std::array<std::string, STAGE2_EXTRA_FEATURES> names = {
      "x_fragment_fraction", "x_fragment_coelution", "x_library_correlation", "x_ms2_abs_ppm", "x_ms1_abs_ppm", "x_ms1_coelution",
      "x_ms1_isotope", "x_ms1_log_area", "x_ms2_log_area", "x_im_centroid_offset", "x_im_coelution", "x_im_central_share"};
    return names;
  }

  std::array<double, STAGE2_EXTRA_FEATURES> stage2ExtraFeatures(const float* tensor, const std::array<double, 6>& library_intensities,
                                                                const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence,
                                                                const Stage2ExtractionOptions& options)
  {
    std::array<double, STAGE2_EXTRA_FEATURES> f{};
    const auto base = stage2Features(tensor, library_intensities, fragment_mz, evidence, options);
    for (std::size_t k = 0; k < 7; ++k) f[k] = base[k + 1];
    const double width = 2.0 * options.rt_half_width / RT_BINS;
    double ms1 = 0.0, ms2 = 0.0;
    // Mobility profile per fragment over the core scan bins.
    std::vector<std::vector<double>> profiles;
    std::array<double, IM_BINS> total{};
    for (std::size_t c = 0; c < FRAGMENTS; ++c)
    {
      if (fragment_mz[c] == 0.0) continue;
      std::vector<double> profile(IM_BINS, 0.0);
      double sum = 0.0;
      for (std::size_t b = 0; b < RT_BINS; ++b)
      {
        const double centre = -options.rt_half_width + (b + 0.5) * width;
        if (std::abs(centre) > options.core_half_width || !(evidence.ms2_mask >> b & 1u)) continue;
        for (std::size_t m = 0; m < IM_BINS; ++m)
        {
          const double v = tensor[(c * RT_BINS + b) * IM_BINS + m];
          profile[m] += v;
          total[m] += v;
          sum += v;
        }
      }
      ms2 += sum;
      if (sum > 0.0) profiles.push_back(std::move(profile));
    }
    for (std::size_t b = 0; b < RT_BINS; ++b)
    {
      const double centre = -options.rt_half_width + (b + 0.5) * width;
      if (std::abs(centre) <= options.core_half_width && (evidence.ms1_mask >> b & 1u)) ms1 += evidence.ms1[b];
    }
    f[7] = std::log1p(ms1);
    f[8] = std::log1p(ms2);
    if (ms2 > 0.0)
    {
      double centroid = 0.0, central = 0.0;
      for (std::size_t m = 0; m < IM_BINS; ++m)
      {
        centroid += total[m] * (m + 0.5);
        if (m >= IM_BINS / 2 - 2 && m < IM_BINS / 2 + 2) central += total[m];
      }
      f[9] = std::abs(centroid / ms2 - IM_BINS / 2.0);
      f[11] = central / ms2;
    }
    else
    {
      f[9] = IM_BINS / 2.0;
    }
    if (profiles.size() >= 2)
    {
      double sum = 0.0;
      std::size_t pairs = 0;
      for (std::size_t i = 0; i < profiles.size(); ++i)
        for (std::size_t j = i + 1; j < profiles.size(); ++j)
        {
          sum += pearson(profiles[i], profiles[j]);
          ++pairs;
        }
      f[10] = sum / pairs;
    }
    return f;
  }

  const std::array<std::string, STAGE2_FAINT_FEATURES>& stage2FaintFeatureNames()
  {
    static const std::array<std::string, STAGE2_FAINT_FEATURES> names = {"p_cos_apex", "p_local_sn"};
    return names;
  }

  std::array<double, STAGE2_FAINT_FEATURES> stage2FaintFeatures(const float* tensor,
                                                                const std::array<double, 6>& library_intensities,
                                                                const std::array<double, 6>& fragment_mz,
                                                                const TraceEvidence& evidence,
                                                                const Stage2ExtractionOptions& options)
  {
    std::array<double, STAGE2_FAINT_FEATURES> f{};
    (void)options;
    std::vector<std::size_t> bins;
    std::vector<std::size_t> channels;
    for (std::size_t b = 0; b < RT_BINS; ++b)
      if (evidence.ms2_mask >> b & 1u) bins.push_back(b);
    for (std::size_t c = 0; c < FRAGMENTS; ++c)
      if (fragment_mz[c] != 0.0) channels.push_back(c);
    if (bins.empty() || channels.empty()) return f;

    std::size_t hint = 0;
    const auto distance = [](std::size_t b) { return b > RT_BINS / 2 ? b - RT_BINS / 2 : RT_BINS / 2 - b; };
    for (std::size_t s = 1; s < bins.size(); ++s)
      if (distance(bins[s]) < distance(bins[hint])) hint = s;

    std::vector<double> consensus(RT_BINS, 0.0);
    std::vector<double> apex_vector(channels.size(), 0.0);
    for (std::size_t c = 0; c < channels.size(); ++c)
    {
      for (const std::size_t b : bins)
      {
        double value = 0.0;
        for (std::size_t m = 0; m < IM_BINS; ++m) value += tensor[(channels[c] * RT_BINS + b) * IM_BINS + m];
        consensus[b] += value;
      }
      for (std::size_t s = hint > 0 ? hint - 1 : 0; s <= std::min(hint + 1, bins.size() - 1); ++s)
      {
        double value = 0.0;
        for (std::size_t m = 0; m < IM_BINS; ++m) value += tensor[(channels[c] * RT_BINS + bins[s]) * IM_BINS + m];
        apex_vector[c] += value;
      }
    }
    double dot = 0.0, norm_x = 0.0, norm_library = 0.0;
    for (std::size_t c = 0; c < channels.size(); ++c)
    {
      const double weight = std::sqrt(std::max(0.0, library_intensities[channels[c]]));
      dot += apex_vector[c] * weight;
      norm_x += apex_vector[c] * apex_vector[c];
      norm_library += weight * weight;
    }
    if (norm_x > 0.0 && norm_library > 0.0) f[static_cast<std::size_t>(Stage2FaintFeature::COS_APEX)] = dot / std::sqrt(norm_x * norm_library);

    std::vector<float> traces(channels.size() * bins.size(), 0.0f);
    for (std::size_t c = 0; c < channels.size(); ++c)
      for (std::size_t s = 0; s < bins.size(); ++s)
        for (std::size_t m = 0; m < IM_BINS; ++m) traces[c * bins.size() + s] += tensor[(channels[c] * RT_BINS + bins[s]) * IM_BINS + m];
    std::vector<double> scratch;
    const PickedPeak peak = pickPeakGroup(traces.data(), channels.size(), bins.size(), bins.size(), 1, hint, PeakPickerOptions{}, scratch);
    if (peak.found)
    {
      std::vector<double> outside;
      for (std::size_t s = 0; s < bins.size(); ++s)
        if (s < peak.left || s > peak.right) outside.push_back(consensus[bins[s]]);
      if (!outside.empty())
      {
        const auto middle = outside.begin() + outside.size() / 2;
        std::nth_element(outside.begin(), middle, outside.end());
        double median = *middle;
        if (outside.size() % 2 == 0)
        {
          const auto lower = std::max_element(outside.begin(), middle);
          median = (*lower + *middle) / 2.0;
        }
        f[static_cast<std::size_t>(Stage2FaintFeature::LOCAL_SN)] = std::log1p(consensus[bins[hint]]) - std::log1p(median);
      }
    }
    return f;
  }

  namespace
  {
    // Median of a copy (average of the two middle values for even length); 0 for empty input.
    double mfMedian(std::vector<double> v)
    {
      if (v.empty()) return 0.0;
      std::sort(v.begin(), v.end());
      const std::size_t n = v.size();
      return n % 2 ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2.0;
    }

    double mfPearson(const std::vector<double>& x, const std::vector<double>& y)
    {
      const double n = static_cast<double>(x.size());
      double mx = 0.0, my = 0.0;
      for (std::size_t i = 0; i < x.size(); ++i) { mx += x[i]; my += y[i]; }
      mx /= n;
      my /= n;
      double sxy = 0.0, sxx = 0.0, syy = 0.0;
      for (std::size_t i = 0; i < x.size(); ++i)
      {
        sxy += (x[i] - mx) * (y[i] - my);
        sxx += (x[i] - mx) * (x[i] - mx);
        syy += (y[i] - my) * (y[i] - my);
      }
      if (sxx <= 0.0 || syy <= 0.0) return 0.0;
      return sxy / std::sqrt(sxx * syy);
    }

    double mfSlog(double x) { return (x < 0.0 ? -1.0 : 1.0) * std::log1p(std::abs(x)); }
  }

  double stage2RobustZ(const std::vector<double>& x, std::size_t centre, std::size_t a, std::size_t b)
  {
    const std::size_t n = x.size();
    std::vector<double> s(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) s[i] = x[i] + (i > 0 ? x[i - 1] : 0.0) + (i + 1 < n ? x[i + 1] : 0.0);
    const auto collect = [&](std::size_t end_left, std::size_t start_right)
    {
      std::vector<double> bg(s.begin(), s.begin() + std::min(end_left, n));
      if (start_right < n) bg.insert(bg.end(), s.begin() + start_right, s.end());
      return bg;
    };
    std::vector<double> bg = collect(a > 1 ? a - 1 : 0, b + 1);
    if (bg.size() < 3) bg = collect(centre > 3 ? centre - 3 : 0, centre + 4);
    if (bg.empty()) return 0.0;
    const double med = mfMedian(bg);
    for (double& v : bg) v = std::abs(v - med);
    const double mad = mfMedian(bg);
    return (s[centre] - med) / (1.4826 * mad + std::sqrt(std::max(med, 0.0)) + 1.0);
  }

  const std::array<std::string, STAGE2_MATCHED_FILTER_FEATURES>& stage2MatchedFilterFeatureNames()
  {
    static const std::array<std::string, STAGE2_MATCHED_FILTER_FEATURES> names = {
      "mt_z_med", "mt_z_trim", "mt_n_z2", "mt_n_z3", "mt_loo_med", "mt_loo_n", "mt_cos", "mt_dom", "mt_llr"};
    return names;
  }

  std::array<double, STAGE2_MATCHED_FILTER_FEATURES> stage2MatchedFilterFeatures(const float* tensor,
                                                                                  const std::array<double, 6>& library_intensities,
                                                                                  const std::array<double, 6>& fragment_mz,
                                                                                  const TraceEvidence& evidence,
                                                                                  const Stage2ExtractionOptions& options)
  {
    using F = Stage2MatchedFilterFeature;
    std::array<double, STAGE2_MATCHED_FILTER_FEATURES> f{};
    std::vector<std::size_t> channels;
    for (std::size_t c = 0; c < FRAGMENTS; ++c)
      if (fragment_mz[c] != 0.0) channels.push_back(c);
    const std::size_t nf = channels.size();
    if (nf < 2) return f;
    const std::size_t ca = RT_BINS / 2;
    const auto at = [&](std::size_t c, std::size_t b, std::size_t m) { return static_cast<double>(tensor[(c * RT_BINS + b) * IM_BINS + m]); };

    // Library weights: sqrt intensities with unit L2 norm, and their squared shares.
    std::vector<double> lw(nf), lshare(nf, 0.0);
    double l2 = 0.0;
    for (std::size_t k = 0; k < nf; ++k)
    {
      lw[k] = std::sqrt(std::max(0.0, library_intensities[channels[k]]));
      l2 += lw[k] * lw[k];
    }
    if (l2 > 0.0)
      for (std::size_t k = 0; k < nf; ++k)
      {
        lw[k] /= std::sqrt(l2);
        lshare[k] = lw[k] * lw[k];
      }

    // Mobility apex of the smoothed core profile and the Gaussian mobility weight.
    std::array<double, IM_BINS> core{};
    double core_sum = 0.0;
    for (const std::size_t c : channels)
      for (std::size_t b = ca - 2; b <= ca + 2; ++b)
        for (std::size_t m = 0; m < IM_BINS; ++m) core[m] += at(c, b, m);
    for (double v : core) core_sum += v;
    std::size_t jb = IM_BINS / 2;
    if (core_sum != 0.0)
    {
      std::array<double, IM_BINS> sm{};
      for (std::size_t m = 0; m < IM_BINS; ++m)
        sm[m] = 0.5 * core[m] + (m > 0 ? 0.25 * core[m - 1] : 0.0) + (m + 1 < IM_BINS ? 0.25 * core[m + 1] : 0.0);
      jb = static_cast<std::size_t>(std::max_element(sm.begin(), sm.end()) - sm.begin());
    }
    const double sd = 0.008 / (2.0 * options.im_half_width / static_cast<double>(IM_BINS));
    std::array<double, IM_BINS> w{};
    for (std::size_t m = 0; m < IM_BINS; ++m)
    {
      const double d = (static_cast<double>(m) - static_cast<double>(jb)) / sd;
      w[m] = std::exp(-0.5 * d * d);
    }
    std::vector<std::vector<double>> X(nf, std::vector<double>(RT_BINS, 0.0));
    for (std::size_t k = 0; k < nf; ++k)
      for (std::size_t b = 0; b < RT_BINS; ++b)
        for (std::size_t m = 0; m < IM_BINS; ++m) X[k][b] += at(channels[k], b, m) * w[m];

    // Peak boundaries as in stage2FaintFeatures (picked around the centre over scanned bins).
    std::size_t a = ca - 1, bnd = ca + 2;
    std::vector<std::size_t> bins;
    for (std::size_t b = 0; b < RT_BINS; ++b)
      if (evidence.ms2_mask >> b & 1u) bins.push_back(b);
    if (!bins.empty())
    {
      std::size_t hint = 0;
      const auto distance = [](std::size_t b) { return b > RT_BINS / 2 ? b - RT_BINS / 2 : RT_BINS / 2 - b; };
      for (std::size_t s = 1; s < bins.size(); ++s)
        if (distance(bins[s]) < distance(bins[hint])) hint = s;
      std::vector<float> traces(nf * bins.size(), 0.0f);
      for (std::size_t c = 0; c < nf; ++c)
        for (std::size_t s = 0; s < bins.size(); ++s)
          for (std::size_t m = 0; m < IM_BINS; ++m) traces[c * bins.size() + s] += tensor[(channels[c] * RT_BINS + bins[s]) * IM_BINS + m];
      std::vector<double> scratch;
      const PickedPeak peak = pickPeakGroup(traces.data(), nf, bins.size(), bins.size(), 1, hint, PeakPickerOptions{}, scratch);
      if (peak.found)
      {
        a = bins[peak.left];
        bnd = bins[peak.right] + 1;
      }
    }

    std::vector<double> z(nf);
    for (std::size_t k = 0; k < nf; ++k) z[k] = stage2RobustZ(X[k], ca, a, bnd);
    f[static_cast<std::size_t>(F::Z_MED)] = mfSlog(mfMedian(z));
    {
      std::vector<double> sorted = z;
      std::sort(sorted.begin(), sorted.end());
      f[static_cast<std::size_t>(F::Z_TRIM)] = mfSlog(std::accumulate(sorted.begin(), sorted.end() - 1, 0.0) / static_cast<double>(nf - 1));
    }
    double llr = 0.0;
    for (std::size_t k = 0; k < nf; ++k)
    {
      if (z[k] >= 2.0) f[static_cast<std::size_t>(F::N_Z2)] += 1.0;
      if (z[k] >= 3.0) f[static_cast<std::size_t>(F::N_Z3)] += 1.0;
      llr += lw[k] * std::log1p(std::max(z[k], 0.0));
    }
    f[static_cast<std::size_t>(F::LLR)] = mfSlog(llr);

    // Leave-one-out correlation over the core RT bins.
    std::vector<double> corrs(nf);
    for (std::size_t k = 0; k < nf; ++k)
    {
      std::vector<double> own, rest;
      for (std::size_t b = ca - 2; b <= ca + 2; ++b)
      {
        own.push_back(X[k][b]);
        double r = 0.0;
        for (std::size_t j = 0; j < nf; ++j)
          if (j != k) r += X[j][b];
        rest.push_back(r);
      }
      corrs[k] = mfPearson(own, rest);
      if (corrs[k] > 0.5) f[static_cast<std::size_t>(F::LOO_N)] += 1.0;
    }
    f[static_cast<std::size_t>(F::LOO_MED)] = mfMedian(corrs);

    // Centre-bin cosine and dominance.
    double dot = 0.0, vv = 0.0, vs = 0.0;
    std::size_t vmax = 0;
    for (std::size_t k = 0; k < nf; ++k)
    {
      const double v = X[k][ca];
      dot += v * lw[k];
      vv += v * v;
      vs += v;
      if (v > X[vmax][ca]) vmax = k;
    }
    if (vv > 0.0 && l2 > 0.0) f[static_cast<std::size_t>(F::COS)] = dot / std::sqrt(vv);
    if (vs != 0.0 && lshare[vmax] > 0.0)
    {
      const double ratio = (X[vmax][ca] / vs) / lshare[vmax];
      if (ratio > 0.0) f[static_cast<std::size_t>(F::DOM)] = std::log(ratio);
    }
    return f;
  }

  const std::array<std::string, STAGE2_MOBILITY_APEX_FEATURES>& stage2MobilityApexFeatureNames()
  {
    static const std::array<std::string, STAGE2_MOBILITY_APEX_FEATURES> names = {
      "m_apex_offset", "m_fragment_sd", "m_apex_cosine", "m_apex_support", "m_apex_contrast"};
    return names;
  }

  std::array<double, STAGE2_MOBILITY_APEX_FEATURES> stage2MobilityApexFeatures(const float* histogram, const std::array<double, 6>& library_intensities,
                                                                              const std::array<double, 6>& fragment_mz)
  {
    constexpr std::size_t B = TraceCollector::MOBILITY_HISTOGRAM_BINS;
    constexpr double W = TraceCollector::MOBILITY_HISTOGRAM_BIN;
    std::array<double, STAGE2_MOBILITY_APEX_FEATURES> f{};
    f[0] = B * W / 2.0;
    f[1] = 0.03;
    std::array<double, B> consensus{};
    double total = 0.0;
    if (histogram != nullptr)
      for (std::size_t c = 0; c < FRAGMENTS; ++c)
        if (fragment_mz[c] != 0.0)
          for (std::size_t b = 0; b < B; ++b)
          {
            consensus[b] += histogram[c * B + b];
            total += histogram[c * B + b];
          }
    if (!(total > 0.0)) return f;
    // Apex: argmax of the [1 2 1]/4-smoothed consensus (zero beyond the ends), first maximum on ties.
    std::size_t apex = 0;
    double best = -1.0;
    for (std::size_t b = 0; b < B; ++b)
    {
      const double s = 0.5 * consensus[b] + 0.25 * (b > 0 ? consensus[b - 1] : 0.0) + 0.25 * (b + 1 < B ? consensus[b + 1] : 0.0);
      if (s > best) { best = s; apex = b; }
    }
    const auto distance = [apex](std::size_t b) { return b > apex ? b - apex : apex - b; };
    // Bin centres relative to the precursor 1/K0: -B*W/2 + (b + 0.5) W.
    const auto centre = [](std::size_t b) { return -static_cast<double>(B) * W / 2.0 + (static_cast<double>(b) + 0.5) * W; };
    f[0] = std::abs(centre(apex));
    std::vector<double> means, weights;
    std::array<double, FRAGMENTS> near{};
    double inside = 0.0, side = 0.0;
    for (std::size_t c = 0; c < FRAGMENTS; ++c)
    {
      if (fragment_mz[c] == 0.0) continue;
      double s = 0.0, m = 0.0;
      for (std::size_t b = 0; b < B; ++b)
      {
        const double v = histogram[c * B + b];
        const std::size_t d = distance(b);
        if (d <= 5) { s += v; m += v * centre(b); }
        if (d <= 2) near[c] += v;
      }
      if (s > 0.0) { means.push_back(m / s); weights.push_back(s); }
    }
    for (std::size_t b = 0; b < B; ++b)
    {
      const std::size_t d = distance(b);
      if (d <= 2) inside += consensus[b];
      else if (d >= 5 && d <= 12) side += consensus[b];
    }
    if (means.size() >= 2)
    {
      double sw = 0.0, mean = 0.0;
      for (std::size_t k = 0; k < means.size(); ++k) { sw += weights[k]; mean += weights[k] * means[k]; }
      mean /= sw;
      double var = 0.0;
      for (std::size_t k = 0; k < means.size(); ++k) var += weights[k] * (means[k] - mean) * (means[k] - mean);
      f[1] = std::sqrt(var / sw);
    }
    // Cosine with sqrt library intensities over the non-padding channels.
    double dot = 0.0, nv = 0.0, nl = 0.0;
    std::size_t support = 0;
    for (std::size_t c = 0; c < FRAGMENTS; ++c)
    {
      if (fragment_mz[c] == 0.0) continue;
      const double l = std::sqrt(std::max(0.0, library_intensities[c]));
      dot += near[c] * l;
      nv += near[c] * near[c];
      nl += l * l;
      if (near[c] > 0.0) ++support;
    }
    f[2] = nv > 0.0 && nl > 0.0 ? dot / std::sqrt(nv * nl) : 0.0;
    f[3] = static_cast<double>(support);
    f[4] = std::log1p(inside) - std::log1p(side);
    return f;
  }

  std::array<double, STAGE2_PEAKGROUP_FEATURES> stage2PeakGroupFeatures(const float* tensor, const std::array<double, 6>& library_intensities,
                                                                        const std::array<double, 6>& fragment_mz,
                                                                        const std::array<std::int8_t, 6>& fragment_charge,
                                                                        const TraceEvidence& evidence, const PeakGroupEvidence& pg,
                                                                        double tolerance_ppm, double floor, PeakGroupSNParts* sn_parts)
  {
    if (sn_parts != nullptr) *sn_parts = PeakGroupSNParts{};
    constexpr int MAX_DELAY = 10;  // OpenSWATH MRMScoring XCORR_MAX_DELAY
    std::array<double, STAGE2_PEAKGROUP_FEATURES> f{};
    f[0] = f[1] = MAX_DELAY;
    const PickedWindow picked = pickWindow(tensor, fragment_mz, evidence);
    const auto& channels = picked.channels;
    const std::size_t m = channels.size();
    if (m == 0) return f;
    // Library weights over the non-padding channels, normalised to sum 1 (uniform if the library has no intensities).
    std::vector<double> w(m);
    for (std::size_t k = 0; k < m; ++k) w[k] = std::max(0.0, library_intensities[channels[k]]);
    const double wsum = std::accumulate(w.begin(), w.end(), 0.0);
    for (auto& v : w) v = wsum > 0.0 ? v / wsum : 1.0 / static_cast<double>(m);

    // Mass deviation of the core-window monoisotope points; it needs no picked peak.
    double massdev = 0.0, massdev_w = 0.0;
    for (std::size_t k = 0; k < m; ++k)
    {
      const std::size_t c = channels[k];
      const double e = pg.mono[c] > 0.0f ? std::min(1.0, std::abs(static_cast<double>(pg.ppm_sum[c]) / pg.mono[c]) / tolerance_ppm) : 1.0;
      massdev += e / static_cast<double>(m);
      massdev_w += w[k] * e;
    }
    f[6] = massdev;
    f[7] = massdev_w;

    if (!picked.peak.found) return f;
    const std::size_t n = picked.n;
    const auto& peak = picked.peak;
    const std::size_t span = peak.right - peak.left + 1;
    std::vector<double> area(m, 0.0);
    double total_all = 0.0;
    std::size_t with_signal = 0;
    for (std::size_t k = 0; k < m; ++k)
    {
      for (std::size_t s = 0; s < n; ++s) total_all += picked.traces[k * n + s];
      for (std::size_t s = peak.left; s <= peak.right; ++s) area[k] += picked.traces[k * n + s];
      if (area[k] > 0.0) ++with_signal;
    }
    if (span < 3 || with_signal < 2) return f;
    const double area_total = std::accumulate(area.begin(), area.end(), 0.0);

    // Cross-correlation (OpenSWATH MRMScoring::initializeXCorrMatrix): standardised traces over the peak, best lag.
    std::vector<std::vector<double>> z(m, std::vector<double>(span));
    for (std::size_t k = 0; k < m; ++k)
    {
      auto& v = z[k];
      for (std::size_t s = 0; s < span; ++s) v[s] = picked.traces[k * n + peak.left + s];
      const double mean = std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(span);
      double var = 0.0;
      for (double x : v) var += (x - mean) * (x - mean);
      double sd = std::sqrt(var / static_cast<double>(span));
      if (mean == 0.0 && sd == 0.0) continue;
      if (sd == 0.0) sd = 1.0;
      for (auto& x : v) x = (x - mean) / sd;
    }
    const int delay = std::min(MAX_DELAY, static_cast<int>(span));
    std::vector<double> lags, shapes;
    double lag_w = 0.0, shape_w = 0.0;
    for (std::size_t i = 0; i < m; ++i)
    {
      for (std::size_t j = i; j < m; ++j)
      {
        int best_lag = -delay;
        double best = -std::numeric_limits<double>::infinity();
        for (int d = -delay; d <= delay; ++d)
        {
          const int start = std::max(0, -d), end = std::min(static_cast<int>(span), static_cast<int>(span) - d);
          double sxy = 0.0;
          for (int s = start; s < end; ++s) sxy += z[i][s] * z[j][s + d];
          sxy /= static_cast<double>(span);
          if (sxy > best) { best = sxy; best_lag = d; }
        }
        const double lag = std::abs(best_lag);
        lags.push_back(lag);
        shapes.push_back(best);
        const double weight = (i == j ? 1.0 : 2.0) * w[i] * w[j];
        lag_w += weight * lag;
        shape_w += weight * best;
      }
    }
    const double lag_mean = std::accumulate(lags.begin(), lags.end(), 0.0) / static_cast<double>(lags.size());
    double lag_var = 0.0;
    for (double l : lags) lag_var += (l - lag_mean) * (l - lag_mean);
    f[0] = lag_mean + (lags.size() > 1 ? std::sqrt(lag_var / static_cast<double>(lags.size() - 1)) : 0.0);
    f[1] = lag_w;
    f[2] = std::accumulate(shapes.begin(), shapes.end(), 0.0) / static_cast<double>(shapes.size());
    f[3] = shape_w;

    // Signal to noise at the consensus apex, against each fragment's median over the scan bins plus the run floor.
    std::vector<double> scratch(n);
    PeakGroupSNParts parts;
    parts.valid = true;
    parts.channels = static_cast<std::uint8_t>(m);
    for (std::size_t k = 0; k < m; ++k)
    {
      for (std::size_t s = 0; s < n; ++s) scratch[s] = picked.traces[k * n + s];
      std::nth_element(scratch.begin(), scratch.begin() + n / 2, scratch.end());
      parts.apex[k] = picked.traces[k * n + peak.apex];
      parts.median[k] = static_cast<float>(scratch[n / 2]);
    }
    f[4] = peakGroupLogSN(parts, floor);
    if (sn_parts != nullptr) *sn_parts = parts;
    f[5] = total_all > 0.0 ? area_total / total_all : 0.0;

    // Fragment isotopes (OpenSWATH DIAScoring::dia_isotope_scores), weighted by the observed area share.
    static thread_local std::unordered_map<long, std::array<double, PeakGroupEvidence::ISOTOPES + 1>> averagine;
    double iso_corr = 0.0, iso_overlap = 0.0;
    for (std::size_t k = 0; k < m; ++k)
    {
      const std::size_t c = channels[k];
      const double a = area[k] / area_total;
      if (a <= 0.0 || pg.mono[c] <= 0.0f) continue;
      const double z_frag = fragment_charge[c] != 0 ? std::abs(static_cast<double>(fragment_charge[c])) : 1.0;
      const long mass = std::lround(fragment_mz[c] * z_frag);
      auto it = averagine.find(mass);
      if (it == averagine.end())
      {
        std::array<double, PeakGroupEvidence::ISOTOPES + 1> theory{};
        // Senko's averagine composition for this mass, coarse isotope pattern (as OpenSWATH dia_isotope_scores).
        EmpiricalFormula composition;
        composition.estimateFromWeightAndComp(static_cast<double>(mass), 4.9384, 7.7583, 1.3577, 1.4773, 0.0417, 0);
        const auto dist = composition.getIsotopeDistribution(CoarseIsotopePatternGenerator(static_cast<int>(theory.size())));
        for (std::size_t i = 0; i < theory.size() && i < dist.size(); ++i) theory[i] = dist[i].getIntensity();
        it = averagine.emplace(mass, theory).first;
      }
      std::vector<double> observed{pg.mono[c]}, theory(it->second.begin(), it->second.end());
      for (float v : pg.isotope[c]) observed.push_back(v);
      iso_corr += a * pearson(observed, theory);
      double larger = 0.0;
      for (float v : pg.before[c]) larger += v > pg.mono[c] ? 1.0 : 0.0;
      iso_overlap += a * larger;
    }
    f[8] = iso_corr;
    f[9] = iso_overlap;
    return f;
  }

  double peakGroupLogSN(const PeakGroupSNParts& parts, double floor)
  {
    if (!parts.valid) return 0.0;
    double sn = 0.0;
    for (std::size_t k = 0; k < parts.channels; ++k) sn += parts.apex[k] / (static_cast<double>(parts.median[k]) + floor);
    return std::log1p(sn / static_cast<double>(parts.channels));
  }

  std::size_t peakGroupFloorStep(std::size_t windows)
  {
    return std::max<std::size_t>(1, (windows + 19999) / 20000);
  }

  void appendPeakGroupFloorValues(const float* t, std::vector<float>& values)
  {
    for (std::size_t c = 0; c < FRAGMENTS; ++c)
    {
      for (std::size_t b = 0; b < RT_BINS; ++b)
      {
        float v = 0.0f;
        for (std::size_t q = 0; q < IM_BINS; ++q) v += t[(c * RT_BINS + b) * IM_BINS + q];
        if (v > 0.0f) values.push_back(v);
      }
    }
  }

  double peakGroupFloorFromValues(std::vector<float>& values)
  {
    if (values.empty()) return 1.0;
    std::nth_element(values.begin(), values.begin() + values.size() / 2, values.end());
    return values[values.size() / 2];
  }

  double peakGroupNoiseFloor(const TraceCollector& traces)
  {
    const std::size_t n = traces.size();
    const std::size_t step = peakGroupFloorStep(n);
    std::vector<float> values;
    std::array<float, TraceCollector::SAMPLE_SIZE> buffer;
    for (std::size_t i = 0; i < n; i += step) appendPeakGroupFloorValues(traces.windowTensor(i, buffer), values);
    return peakGroupFloorFromValues(values);
  }

  std::array<double, STAGE2_SEARCH_FEATURES> stage2SearchFeatures(const SearchPeakFeatures& peak)
  {
    return {peak.transition_support, std::log1p(peak.matched_scans), std::log1p(peak.peak_width_seconds), peak.abs_rt_error_seconds};
  }

  std::unordered_map<std::string, SearchPeakFeatures> readSearchPeakFeatures(const std::string& path)
  {
    std::ifstream input(path);
    if (!input.good()) throw std::invalid_argument("Could not open search features: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("Search feature table is empty: " + path);
    const auto split = [](std::string text)
    {
      if (!text.empty() && text.back() == '\r') text.pop_back();
      std::vector<std::string> out;
      std::size_t start = 0;
      for (std::size_t tab; (tab = text.find('\t', start)) != std::string::npos; start = tab + 1) out.push_back(text.substr(start, tab - start));
      out.push_back(text.substr(start));
      return out;
    };
    const auto header = split(line);
    const std::array<std::string, 8> names = {"candidate_row", "precursor_id", "rt_apex", "peak_width_seconds", "transition_support",
                                             "matched_scans", "has_library_rt", "apex_rt_error_seconds"};
    std::array<std::size_t, 8> at{};
    for (std::size_t k = 0; k < names.size(); ++k)
    {
      const auto found = std::find(header.begin(), header.end(), names[k]);
      if (found == header.end() || std::count(header.begin(), header.end(), names[k]) != 1)
      {
        throw std::invalid_argument("Search feature table needs exactly one column " + names[k]);
      }
      at[k] = static_cast<std::size_t>(found - header.begin());
    }
    const auto number = [](const std::string& text)
    {
      std::size_t used = 0;
      const double value = std::stod(text, &used);
      if (used != text.size() || !std::isfinite(value)) throw std::invalid_argument("Search feature value is not a finite number: " + text);
      return value;
    };
    std::unordered_map<std::string, SearchPeakFeatures> out;
    while (std::getline(input, line))
    {
      const auto values = split(line);
      if (values.size() != header.size()) throw std::invalid_argument("Search feature row has a different column count than the header.");
      SearchPeakFeatures peak;
      peak.precursor_id = values[at[1]];
      peak.rt_apex = number(values[at[2]]);
      peak.peak_width_seconds = number(values[at[3]]);
      peak.transition_support = number(values[at[4]]);
      peak.matched_scans = number(values[at[5]]);
      peak.abs_rt_error_seconds = number(values[at[6]]) != 0.0 ? std::abs(number(values[at[7]])) : 0.0;
      if (peak.peak_width_seconds < 0.0 || peak.transition_support < 0.0 || peak.matched_scans < 0.0)
      {
        throw std::invalid_argument("Search feature counts and widths must be non-negative.");
      }
      if (!out.emplace(values[at[0]], std::move(peak)).second) throw std::invalid_argument("Duplicate candidate_row in search features: " + values[at[0]]);
    }
    return out;
  }

  std::vector<double> crossfitStage2LDA(const std::vector<Stage2Row>& rows, int folds, const Stage2Options& options,
                                        std::vector<Stage2FoldReport>& report)
  {
    if (rows.empty() || rows.front().features.empty())
    {
      throw std::invalid_argument("Second-stage LDA needs rows with at least one feature.");
    }
    const std::size_t P = rows.front().features.size();
    if (folds < 2 || options.iterations < 1 || !(options.initial_fdr > 0.0) || !(options.iteration_fdr > 0.0) ||
        !(options.ridge >= 0.0) || options.min_positives < 2)
    {
      throw std::invalid_argument("Second-stage LDA needs at least two folds, one iteration, positive FDR cuts and a non-negative ridge.");
    }
    for (const auto& row : rows)
    {
      if (row.fold < 0 || row.fold >= folds || row.precursor_id.empty() || row.features.size() != P)
      {
        throw std::invalid_argument("Second-stage rows need a precursor id, a fold in 0..folds-1 and the same feature count.");
      }
      for (double v : row.features)
      {
        if (!std::isfinite(v)) throw std::invalid_argument("Second-stage features must be finite.");
      }
    }
    std::vector<double> out(rows.size(), 0.0);
    report.clear();
    for (int fold = 0; fold < folds; ++fold)
    {
      std::vector<std::size_t> train, test;
      for (std::size_t i = 0; i < rows.size(); ++i) (rows[i].fold == fold ? test : train).push_back(i);
      if (train.empty() || test.empty()) throw std::invalid_argument("Every second-stage fold needs training and held-out rows.");

      // Standardise with the training rows' mean and SD (a constant feature keeps SD 1).
      std::vector<double> mean(P, 0.0), sd(P, 0.0);
      for (std::size_t i : train) for (std::size_t p = 0; p < P; ++p) mean[p] += rows[i].features[p];
      for (std::size_t p = 0; p < P; ++p) mean[p] /= train.size();
      for (std::size_t i : train)
        for (std::size_t p = 0; p < P; ++p) sd[p] += (rows[i].features[p] - mean[p]) * (rows[i].features[p] - mean[p]);
      for (std::size_t p = 0; p < P; ++p)
      {
        sd[p] = std::sqrt(sd[p] / train.size());
        if (!(sd[p] > 0.0)) sd[p] = 1.0;
      }
      const auto z = [&](std::size_t i, std::size_t p) { return (rows[i].features[p] - mean[p]) / sd[p]; };

      std::vector<double> score(rows.size(), 0.0);
      for (std::size_t i : train) score[i] = rows[i].features[0];
      Stage2FoldReport fold_report;
      fold_report.fold = fold;
      std::vector<double> weights(P, 0.0);
      for (int iteration = 0; iteration < options.iterations; ++iteration)
      {
        std::vector<std::size_t> positives, negatives;
        selectTraining(rows, winners(rows, train, score), score, iteration == 0 ? options.initial_fdr : options.iteration_fdr,
                       positives, negatives);
        if (positives.size() < options.min_positives || negatives.size() < options.min_positives)
        {
          throw std::invalid_argument("Second-stage fold " + std::to_string(fold) + " iteration " + std::to_string(iteration) +
                                      " has " + std::to_string(positives.size()) + " positives and " +
                                      std::to_string(negatives.size()) + " negatives; too few to fit.");
        }
        // Ridge LDA: w = (S_within + ridge I)^-1 (mu_pos - mu_neg) on standardised features.
        std::vector<double> mp(P, 0.0), mn(P, 0.0);
        for (std::size_t i : positives) for (std::size_t p = 0; p < P; ++p) mp[p] += z(i, p);
        for (std::size_t i : negatives) for (std::size_t p = 0; p < P; ++p) mn[p] += z(i, p);
        for (std::size_t p = 0; p < P; ++p)
        {
          mp[p] /= positives.size();
          mn[p] /= negatives.size();
        }
        std::vector<double> s(P * P, 0.0);
        const auto add = [&](const std::vector<std::size_t>& set, const std::vector<double>& mu)
        {
          std::vector<double> d(P);
          for (std::size_t i : set)
          {
            for (std::size_t p = 0; p < P; ++p) d[p] = z(i, p) - mu[p];
            for (std::size_t p = 0; p < P; ++p) for (std::size_t q = 0; q < P; ++q) s[p * P + q] += d[p] * d[q];
          }
        };
        add(positives, mp);
        add(negatives, mn);
        const double n = static_cast<double>(positives.size() + negatives.size() - 2);
        for (auto& v : s) v /= n;
        for (std::size_t p = 0; p < P; ++p) s[p * P + p] += options.ridge;
        std::vector<double> diff(P);
        for (std::size_t p = 0; p < P; ++p) diff[p] = mp[p] - mn[p];
        weights = solveSPD(s, diff, P);
        for (std::size_t i : train)
        {
          double v = 0.0;
          for (std::size_t p = 0; p < P; ++p) v += weights[p] * z(i, p);
          score[i] = v;
        }
        fold_report.positives = positives.size();
        fold_report.negatives = negatives.size();
        fold_report.iterations = iteration + 1;
      }
      fold_report.weights = weights;
      for (std::size_t i : test)
      {
        double v = 0.0;
        for (std::size_t p = 0; p < P; ++p) v += weights[p] * z(i, p);
        out[i] = v;
      }
      report.push_back(fold_report);
    }
    return out;
  }
}

namespace OpenMS::OpenDIA
{
  std::vector<float> stage2SiblingFeatures(const std::vector<std::uint32_t>& group, const std::vector<std::uint32_t>& precursor,
                                           const std::vector<double>& rt, const std::vector<float>& score,
                                           const std::vector<std::uint32_t>& group_precursors, float not_kept, double rt_tolerance)
  {
    const std::size_t n = group.size();
    constexpr std::size_t F = STAGE2_SIBLING_FEATURES;
    float smin = std::numeric_limits<float>::max();
    bool any = false;
    for (float s : score)
      if (s > not_kept) { smin = std::min(smin, s); any = true; }
    if (!any) smin = 0.0f;
    // Windows sorted by (group, rt, row): each group is one contiguous run.
    std::vector<std::uint32_t> order(n);
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b)
              { return group[a] != group[b] ? group[a] < group[b] : (rt[a] != rt[b] ? rt[a] < rt[b] : a < b); });
    std::vector<float> out(n * F, 0.0f);
    std::vector<std::uint32_t> seen;
    std::size_t lo = 0;
    for (std::size_t k = 0; k < n; ++k)
    {
      const std::uint32_t i = order[k];
      while (group[order[lo]] != group[i] || rt[i] - rt[order[lo]] > rt_tolerance) ++lo;
      float best = smin;
      double best_d = rt_tolerance;
      bool found = false;
      seen.clear();
      for (std::size_t m = lo; m < n && group[order[m]] == group[i] && rt[order[m]] - rt[i] <= rt_tolerance; ++m)
      {
        const std::uint32_t j = order[m];
        if (precursor[j] == precursor[i]) continue;
        if (std::find(seen.begin(), seen.end(), precursor[j]) == seen.end()) seen.push_back(precursor[j]);
        const float s = score[j] > not_kept ? score[j] : smin;
        const double d = std::abs(rt[j] - rt[i]);
        if (!found || s > best || (s == best && d < best_d)) { best = s; best_d = d; found = true; }
      }
      float* f = out.data() + i * F;
      f[0] = static_cast<float>(group_precursors[group[i]] > 0 ? group_precursors[group[i]] - 1 : 0);
      f[1] = static_cast<float>(seen.size());
      f[2] = best;
      f[3] = static_cast<float>(best_d);
    }
    return out;
  }
}
