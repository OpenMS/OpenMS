// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIAStage2.h"
#include "OpenDIAPeakPicker.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

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

  const std::array<std::string, STAGE2_SEARCH_FEATURES>& stage2SearchFeatureNames()
  {
    static const std::array<std::string, STAGE2_SEARCH_FEATURES> names = {
      "transition_support", "log_matched_scans", "log_peak_width", "abs_library_rt_error"};
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

    // Peak boundaries picked around the centre over scanned bins.
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

  std::array<double, STAGE2_SEARCH_FEATURES> stage2SearchFeatures(const SearchPeakFeatures& peak)
  {
    return {peak.transition_support, std::log1p(peak.matched_scans), std::log1p(peak.peak_width_seconds), peak.abs_rt_error_seconds};
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
