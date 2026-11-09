// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIAStage2.h"
#include "OpenDIAPeakPicker.h"

#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>

#include <algorithm>
#include <cmath>
#include <fstream>
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
    Stage2FeatureSet set{false, false, false};
    const std::string spec = name == "extraction_search" ? std::string("extraction+search") : name;
    std::size_t start = 0;
    while (true)
    {
      const std::size_t plus = spec.find('+', start);
      const std::string group = spec.substr(start, plus == std::string::npos ? std::string::npos : plus - start);
      bool* flag = group == "extraction" ? &set.extraction : group == "search" ? &set.search : group == "isotope" ? &set.isotope
                   : group == "boundary" ? &set.boundary : group == "byion" ? &set.byion : group == "shape" ? &set.shape : nullptr;
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
    if (set.isotope)
    {
      names.insert(names.end(), stage2IsotopeFeatureNames().begin(), stage2IsotopeFeatureNames().end());
    }
    if (set.boundary)
    {
      names.insert(names.end(), stage2BoundaryFeatureNames().begin(), stage2BoundaryFeatureNames().end());
    }
    if (set.byion)
    {
      names.insert(names.end(), stage2ByIonFeatureNames().begin(), stage2ByIonFeatureNames().end());
    }
    if (set.shape)
    {
      names.insert(names.end(), stage2ShapeFeatureNames().begin(), stage2ShapeFeatureNames().end());
    }
    return names;
  }

  const std::array<std::string, STAGE2_SHAPE_FEATURES>& stage2ShapeFeatureNames()
  {
    static const std::array<std::string, STAGE2_SHAPE_FEATURES> names = {
      "shape_xcorr_coelution", "shape_xcorr_shape", "shape_xcorr_coelution_w", "shape_xcorr_shape_w", "shape_xcorr_shape_min",
      "shape_library_dotproduct", "shape_library_manhattan", "shape_library_rmsd", "shape_library_spectral_angle",
      "shape_log_sn", "shape_intensity_share", "shape_signal_fragments", "shape_min_fragment_share"};
    return names;
  }

  std::array<double, STAGE2_SHAPE_FEATURES> stage2ShapeFeatures(const float* tensor, const std::array<double, 6>& library_intensities,
                                                                const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence)
  {
    std::array<double, STAGE2_SHAPE_FEATURES> f{};
    const PickedWindow picked = pickWindow(tensor, fragment_mz, evidence);
    const std::size_t c = picked.channels.size();
    if (!picked.peak.found || c < 2) return f;
    const std::size_t n = picked.n, left = picked.peak.left, right = picked.peak.right, span = right - left + 1;
    const auto at = [&](std::size_t k, std::size_t s) { return static_cast<double>(picked.traces[k * n + s]); };
    // Library weights and per-fragment window statistics.
    std::vector<double> lib(c), libw(c), area(c, 0.0), median(c), peak_max(c, 0.0);
    double lib_sum = 0.0, window_total = 0.0;
    for (std::size_t k = 0; k < c; ++k)
    {
      lib[k] = std::sqrt(std::max(0.0, library_intensities[picked.channels[k]]));
      lib_sum += lib[k];
      std::vector<double> all(n);
      for (std::size_t s = 0; s < n; ++s)
      {
        all[s] = at(k, s);
        window_total += all[s];
      }
      for (std::size_t s = left; s <= right; ++s)
      {
        area[k] += all[s];
        peak_max[k] = std::max(peak_max[k], all[s]);
      }
      std::nth_element(all.begin(), all.begin() + n / 2, all.end());
      median[k] = all[n / 2];
    }
    for (std::size_t k = 0; k < c; ++k) libw[k] = lib_sum > 0.0 ? lib[k] / lib_sum : 1.0 / c;
    // Lagged cross-correlation over all fragment pairs inside the peak.
    if (span >= 3)
    {
      std::vector<std::vector<double>> centred(c, std::vector<double>(span));
      std::vector<double> norm(c, 0.0);
      for (std::size_t k = 0; k < c; ++k)
      {
        const double mean = area[k] / span;
        for (std::size_t s = 0; s < span; ++s)
        {
          centred[k][s] = at(k, left + s) - mean;
          norm[k] += centred[k][s] * centred[k][s];
        }
        norm[k] = std::sqrt(norm[k]);
      }
      const int max_lag = 3;
      double lag_sum = 0.0, shape_sum = 0.0, lag_w = 0.0, shape_w = 0.0, weight = 0.0, shape_min = 2.0;
      std::size_t pairs = 0;
      for (std::size_t i = 0; i < c; ++i)
      {
        for (std::size_t j = i + 1; j < c; ++j)
        {
          double best = 0.0, best_lag = max_lag;
          if (norm[i] > 0.0 && norm[j] > 0.0)
          {
            best = -2.0;
            for (int lag = -max_lag; lag <= max_lag; ++lag)
            {
              double t = 0.0;
              for (std::size_t s = 0; s < span; ++s)
              {
                const long u = static_cast<long>(s) + lag;
                if (u >= 0 && u < static_cast<long>(span)) t += centred[i][s] * centred[j][u];
              }
              const double v = t / (norm[i] * norm[j]);
              if (v > best)
              {
                best = v;
                best_lag = std::abs(lag);
              }
            }
          }
          const double w = libw[i] * libw[j];
          lag_sum += best_lag;
          shape_sum += best;
          lag_w += w * best_lag;
          shape_w += w * best;
          weight += w;
          shape_min = std::min(shape_min, best);
          ++pairs;
        }
      }
      f[0] = lag_sum / pairs;
      f[1] = shape_sum / pairs;
      f[2] = weight > 0.0 ? lag_w / weight : f[0];
      f[3] = weight > 0.0 ? shape_w / weight : f[1];
      f[4] = shape_min;
    }
    // Library agreement of sqrt areas.
    std::vector<double> root(c);
    double root_sum = 0.0, root_norm = 0.0, lib_norm = 0.0, dot = 0.0, total_area = 0.0;
    for (std::size_t k = 0; k < c; ++k)
    {
      root[k] = std::sqrt(std::max(0.0, area[k]));
      root_sum += root[k];
      root_norm += root[k] * root[k];
      lib_norm += lib[k] * lib[k];
      dot += root[k] * lib[k];
      total_area += std::max(0.0, area[k]);
    }
    if (root_norm > 0.0 && lib_norm > 0.0)
    {
      f[5] = dot / std::sqrt(root_norm * lib_norm);
      double manhattan = 0.0, squares = 0.0;
      for (std::size_t k = 0; k < c; ++k)
      {
        const double d = root[k] / root_sum - libw[k];
        manhattan += std::abs(d);
        squares += d * d;
      }
      f[6] = manhattan;
      f[7] = std::sqrt(squares / c);
      f[8] = 1.0 - 2.0 * std::acos(std::clamp(f[5], -1.0, 1.0)) / 3.14159265358979323846;
    }
    else
    {
      f[6] = 2.0;
      f[7] = std::sqrt(1.0 / c);
    }
    // Signal: apex S/N against the window median, peak share of the window signal, fragments with signal, weakest fragment.
    const std::size_t apex = picked.peak.apex;
    double log_sn = 0.0;
    std::size_t signal = 0;
    double min_area = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < c; ++k)
    {
      const double noise = std::max(median[k], 1e-3 * std::max(peak_max[k], 1.0));
      log_sn += std::log1p(at(k, apex) / noise);
      if (peak_max[k] > 3.0 * noise) ++signal;
      min_area = std::min(min_area, std::max(0.0, area[k]));
    }
    f[9] = log_sn / c;
    f[10] = window_total > 0.0 ? total_area / window_total : 0.0;
    f[11] = static_cast<double>(signal);
    f[12] = total_area > 0.0 ? min_area / total_area : 0.0;
    return f;
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

  const std::array<std::string, STAGE2_ISOTOPE_FEATURES>& stage2IsotopeFeatureNames()
  {
    static const std::array<std::string, STAGE2_ISOTOPE_FEATURES> names = {"ms1_isotope_correlation", "ms1_isotope_overlap", "ms1_isotope_coelution"};
    return names;
  }

  std::array<double, STAGE2_ISOTOPE_FEATURES> stage2IsotopeFeatures(const TraceEvidence& evidence, double precursor_mz, int charge,
                                                                    const Stage2ExtractionOptions& options)
  {
    std::array<double, STAGE2_ISOTOPE_FEATURES> f{};
    if (charge <= 0 || !(precursor_mz > 0.0)) return f;
    // Averagine pattern at the precursor mass (OpenSWATH uses the same model for its MS1 isotope scores).
    const auto& averagine = DIAHelpers::averagineIsotopeDistribution(precursor_mz * charge, 4);
    std::vector<double> observed, expected;
    for (std::size_t k = 0; k < 4; ++k)
    {
      observed.push_back(evidence.ms1_isotope_core[k + 1]);
      expected.push_back(k < averagine.size() ? averagine[k].getIntensity() : 0.0);
    }
    f[0] = observed[0] > 0.0 ? pearson(observed, expected) : 0.0;
    const double mono = evidence.ms1_isotope_core[1];
    f[1] = mono > 0.0 ? std::min(10.0, evidence.ms1_isotope_core[0] / mono) : (evidence.ms1_isotope_core[0] > 0.0 ? 10.0 : 0.0);
    // M+1 and M+2 traces against M+0 over the core bins that had an MS1 scan.
    const double width = 2.0 * options.rt_half_width / RT_BINS;
    std::vector<double> m0, m1, m2;
    for (std::size_t b = 0; b < RT_BINS; ++b)
    {
      const double centre = -options.rt_half_width + (b + 0.5) * width;
      if (std::abs(centre) > options.core_half_width || !(evidence.ms1_mask >> b & 1u)) continue;
      m0.push_back(evidence.ms1[b]);
      m1.push_back(evidence.ms1_isotope_traces[0][b]);
      m2.push_back(evidence.ms1_isotope_traces[1][b]);
    }
    if (m0.size() >= 3) f[2] = 0.5 * (pearson(m0, m1) + pearson(m0, m2));
    return f;
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
