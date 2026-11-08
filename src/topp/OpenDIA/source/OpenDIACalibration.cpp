// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIACalibration.h>

#include <OpenMS/IONMOBILITY/IMTypes.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <unordered_map>

namespace OpenMS::OpenDIA
{
  namespace
  {
    constexpr double PI = 3.14159265358979323846;
    constexpr std::size_t CORRELATED_FRAGMENTS = 4;

    /// About ten anchors per RT-map bin, between 3 and 20 bins.
    std::size_t mapBins_(std::size_t anchors)
    {
      return std::clamp<std::size_t>(anchors / 10, 3, 20);
    }

    double quantile_(std::vector<double> values, double q)
    {
      if (values.empty())
      {
        return 0.0;
      }
      std::sort(values.begin(), values.end());
      return values[static_cast<std::size_t>(std::clamp(q, 0.0, 1.0) * static_cast<double>(values.size() - 1))];
    }

    double median_(std::vector<double> values)
    {
      if (values.empty())
      {
        return 0.0;
      }
      const std::size_t mid = values.size() / 2;
      std::nth_element(values.begin(), values.begin() + mid, values.end());
      const double upper = values[mid];
      if (values.size() % 2 == 1)
      {
        return upper;
      }
      return (*std::max_element(values.begin(), values.begin() + mid) + upper) / 2.0;
    }

    /// 1.4826 x median absolute deviation about @p centre.
    double robustSigma_(const std::vector<double>& values, double centre)
    {
      std::vector<double> deviations;
      deviations.reserve(values.size());
      for (double value : values)
      {
        deviations.push_back(std::fabs(value - centre));
      }
      return 1.4826 * median_(std::move(deviations));
    }

    double pearson_(const std::vector<double>& x, const std::vector<double>& y)
    {
      const double n = static_cast<double>(x.size());
      const double mx = std::accumulate(x.begin(), x.end(), 0.0) / n;
      const double my = std::accumulate(y.begin(), y.end(), 0.0) / n;
      double sxy = 0.0, sxx = 0.0, syy = 0.0;
      for (std::size_t i = 0; i < x.size(); ++i)
      {
        sxy += (x[i] - mx) * (y[i] - my);
        sxx += (x[i] - mx) * (x[i] - mx);
        syy += (y[i] - my) * (y[i] - my);
      }
      return sxx > 0.0 && syy > 0.0 ? sxy / std::sqrt(sxx * syy) : std::numeric_limits<double>::quiet_NaN();
    }

    /// Piecewise-linear interpolation through knots (constant outside).
    double interpolate_(const std::vector<double>& xs, const std::vector<double>& ys, double x)
    {
      if (x <= xs.front())
      {
        return ys.front();
      }
      if (x >= xs.back())
      {
        return ys.back();
      }
      const std::size_t i = std::upper_bound(xs.begin(), xs.end(), x) - xs.begin();
      const double t = xs[i] == xs[i - 1] ? 0.0 : (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
      return ys[i - 1] + t * (ys[i] - ys[i - 1]);
    }

    /// Weighted pool-adjacent-violators: nondecreasing fit of @p y with weights @p w.
    std::vector<double> isotonic_(const std::vector<double>& y, const std::vector<double>& w)
    {
      std::vector<double> value, weight;
      std::vector<std::size_t> size;
      for (std::size_t i = 0; i < y.size(); ++i)
      {
        value.push_back(y[i]);
        weight.push_back(std::max(w[i], 1e-12));
        size.push_back(1);
        while (value.size() > 1 && value[value.size() - 2] > value.back())
        {
          const double total = weight[weight.size() - 2] + weight.back();
          value[value.size() - 2] = (value[value.size() - 2] * weight[weight.size() - 2] + value.back() * weight.back()) / total;
          weight[weight.size() - 2] = total;
          size[size.size() - 2] += size.back();
          value.pop_back();
          weight.pop_back();
          size.pop_back();
        }
      }
      std::vector<double> out;
      for (std::size_t block = 0; block < value.size(); ++block)
      {
        out.insert(out.end(), size[block], value[block]);
      }
      return out;
    }

    /// Paired target-decoy competition; q-values for targets (decoys and pair losers get 1).
    std::vector<double> pairedQValues_(const std::vector<double>& score, const std::vector<bool>& is_decoy,
                                       const std::vector<std::size_t>& partner, const std::vector<bool>& present)
    {
      std::vector<std::size_t> competing;
      for (std::size_t i = 0; i < score.size(); ++i)
      {
        if (!present[i])
        {
          continue;
        }
        const std::size_t p = partner[i];
        if (!present[p])
        {
          competing.push_back(i);
          continue;
        }
        if (is_decoy[i])
        {
          continue; // Decided from the target side.
        }
        // Exact ties go to the decoy.
        competing.push_back(score[i] > score[p] ? i : p);
      }
      std::sort(competing.begin(), competing.end(), [&](std::size_t a, std::size_t b)
      {
        return score[a] != score[b] ? score[a] > score[b] : a < b;
      });
      std::vector<double> q(score.size(), 1.0);
      std::vector<double> fdr(competing.size(), 1.0);
      std::size_t targets = 0, decoys = 0;
      for (std::size_t begin = 0; begin < competing.size();)
      {
        std::size_t end = begin;
        while (end < competing.size() && score[competing[end]] == score[competing[begin]])
        {
          (is_decoy[competing[end]] ? decoys : targets) += 1;
          ++end;
        }
        const double value = targets == 0 ? 1.0 : std::min(1.0, (decoys + 1.0) / static_cast<double>(targets));
        std::fill(fdr.begin() + static_cast<std::ptrdiff_t>(begin), fdr.begin() + static_cast<std::ptrdiff_t>(end), value);
        begin = end;
      }
      double minimum = 1.0;
      for (std::size_t k = competing.size(); k-- > 0;)
      {
        minimum = std::min(minimum, fdr[k]);
        if (!is_decoy[competing[k]])
        {
          q[competing[k]] = minimum;
        }
      }
      return q;
    }
  }

  CalibrationCollector::CalibrationCollector(const OpenSwath::LightTargetedExperiment& library,
                                             const std::vector<std::pair<std::string, std::string>>& traml_ids,
                                             CalibrationOptions options, std::optional<CalibrationWindow> window) :
    options_(options), window_(std::move(window))
  {
    if (options_.panel_targets < 50 || !(options_.search_ppm > 0.0) || !(options_.search_im > 0.0) ||
        options_.candidates == 0 || options_.seed_half_window == 0 || options_.em_knots < 3 ||
        !(options_.anchor_fdr > 0.0 && options_.anchor_fdr < 1.0))
    {
      throw std::invalid_argument("Invalid calibration options");
    }
    // Pair DECOY_<X> with <X> by TRAML_ID.
    std::unordered_map<std::string, std::string> name_of(traml_ids.begin(), traml_ids.end());
    std::unordered_map<std::string, std::string> id_of_name;
    for (const auto& [id, name] : traml_ids)
    {
      id_of_name.emplace(name, id);
    }
    std::unordered_map<std::string, std::vector<const OpenSwath::LightTransition*>> transitions;
    for (const auto& transition : library.getTransitions())
    {
      transitions[transition.peptide_ref].push_back(&transition);
    }
    std::unordered_map<std::string, const OpenSwath::LightCompound*> compounds;
    for (const auto& compound : library.getCompounds())
    {
      compounds.emplace(compound.id, &compound);
    }
    const auto make = [&](const std::string& id, Member& member)
    {
      const auto found = transitions.find(id);
      const auto compound = compounds.find(id);
      if (found == transitions.end() || compound == compounds.end() || !std::isfinite(compound->second->rt))
      {
        return false;
      }
      auto fragments = found->second;
      std::stable_sort(fragments.begin(), fragments.end(), [](const auto* a, const auto* b)
      {
        return a->library_intensity != b->library_intensity ? a->library_intensity > b->library_intensity : a->product_mz < b->product_mz;
      });
      if (fragments.size() < options_.strict_min_fragments)
      {
        return false;
      }
      member.id = id;
      member.precursor_mz = fragments.front()->precursor_mz;
      member.library_rt = compound->second->rt;
      member.library_im = std::isfinite(compound->second->drift_time) && compound->second->drift_time > 0.0 ? compound->second->drift_time : 0.0;
      member.is_decoy = fragments.front()->getDecoy();
      member.fragment_count = std::min(FRAGMENTS, fragments.size());
      for (std::size_t f = 0; f < member.fragment_count; ++f)
      {
        member.fragment_mz[f] = fragments[f]->product_mz;
        member.library_intensity[f] = fragments[f]->library_intensity;
      }
      return true;
    };
    // Eligible pairs, sorted by target library RT (ties: target ID).
    std::vector<std::pair<Member, Member>> eligible;
    for (const auto& compound : library.getCompounds())
    {
      const auto name = name_of.find(compound.id);
      if (name == name_of.end() || name->second.rfind("DECOY_", 0) == 0)
      {
        continue;
      }
      const auto decoy = id_of_name.find("DECOY_" + name->second);
      Member target, decoy_member;
      if (decoy == id_of_name.end() || !make(compound.id, target) || target.is_decoy || !make(decoy->second, decoy_member) ||
          !decoy_member.is_decoy)
      {
        continue;
      }
      eligible.emplace_back(std::move(target), std::move(decoy_member));
    }
    if (eligible.size() < 50)
    {
      throw std::invalid_argument("Calibration needs at least 50 target-decoy pairs with library RT and five fragments");
    }
    std::sort(eligible.begin(), eligible.end(), [](const auto& a, const auto& b)
    {
      return a.first.library_rt != b.first.library_rt ? a.first.library_rt < b.first.library_rt : a.first.id < b.first.id;
    });
    // Even stride over library RT.
    const std::size_t count = std::min(options_.panel_targets, eligible.size());
    for (std::size_t i = 0; i < count; ++i)
    {
      auto& pair = eligible[i * eligible.size() / count];
      const std::size_t t = panel_.size();
      pair.first.partner = t + 1;
      pair.second.partner = t;
      panel_.push_back(std::move(pair.first));
      panel_.push_back(std::move(pair.second));
    }
    by_mz_.resize(panel_.size());
    std::iota(by_mz_.begin(), by_mz_.end(), 0);
    std::sort(by_mz_.begin(), by_mz_.end(), [this](std::size_t a, std::size_t b) { return panel_[a].precursor_mz < panel_[b].precursor_mz; });
    if (window_)
    {
      // Windowed pass: keep only spectra inside the stage-3 window around the stage-1 map, plus a margin that
      // leaves room for the co-elution half window and the local-maximum test at the window edges.
      const double margin = 30.0 + 2.0 * window_->peak_fwhm_s;
      for (auto& member : panel_)
      {
        const double centre = window_->map(member.is_decoy ? panel_[member.partner].library_rt : member.library_rt);
        member.window_low = centre - window_->half_window_s - margin;
        member.window_high = centre + window_->half_window_s + margin;
      }
    }
  }

  void CalibrationCollector::consumeSpectrum(MSSpectrum& spectrum)
  {
    if (spectrum.getMSLevel() != 2 || spectrum.getPrecursors().empty() || !std::isfinite(spectrum.getRT()))
    {
      return;
    }
    const double rt = spectrum.getRT();
    if (!saw_spectrum_)
    {
      first_rt_ = rt;
      saw_spectrum_ = true;
    }
    first_rt_ = std::min(first_rt_, rt);
    last_rt_ = std::max(last_rt_, rt);
    if (!spectrum.isSorted())
    {
      spectrum.sortByPosition();
    }
    const DataArrays::FloatDataArray* mobilities = nullptr;
    if (spectrum.containsIMData())
    {
      const auto [index, unit] = spectrum.getIMData();
      if (spectrum.getFloatDataArrays()[index].size() == spectrum.size() &&
          (unit == DriftTimeUnit::VSSC || (ion_mobility_array_as_vssc_ && unit != DriftTimeUnit::CCS)))
      {
        mobilities = &spectrum.getFloatDataArrays()[index];
      }
    }
    const double scalar_im = spectrum.getDriftTime();
    const bool scalar_valid = spectrum.getDriftTimeUnit() == DriftTimeUnit::VSSC && std::isfinite(scalar_im) && scalar_im > 0.0;
    const bool has_im = mobilities != nullptr || scalar_valid;
    saw_im_ = saw_im_ || has_im;

    for (const auto& precursor : spectrum.getPrecursors())
    {
      const double lower = precursor.getMZ() - precursor.getIsolationWindowLowerOffset();
      const double upper = precursor.getMZ() + precursor.getIsolationWindowUpperOffset();
      if (!(upper > lower))
      {
        continue;
      }
      auto it = std::lower_bound(by_mz_.begin(), by_mz_.end(), lower, [this](std::size_t i, double mz) { return panel_[i].precursor_mz < mz; });
      for (; it != by_mz_.end() && panel_[*it].precursor_mz <= upper; ++it)
      {
        Member& member = panel_[*it];
        if (rt < member.window_low || rt > member.window_high)
        {
          continue; // Outside the windowed pass's collection range (never in the default whole-run pass).
        }
        if (!member.rt.empty() && rt < member.rt.back())
        {
          throw std::runtime_error("Calibration requires MS2 spectra in nondecreasing RT order");
        }
        // Spectra sharing one RT (diaPASEF windows of one frame) merge into one entry by the per-fragment maximum.
        const bool merge = !member.rt.empty() && member.rt.back() == rt;
        std::array<float, FRAGMENTS> intensity{}, ppm{}, im{};
        for (std::size_t f = 0; f < member.fragment_count; ++f)
        {
          const double target = member.fragment_mz[f];
          const double tolerance = target * options_.search_ppm * 1e-6;
          // Trace intensity is the sum over all points inside the m/z and 1/K0 tolerance (diaPASEF spreads one
          // fragment over many mobility scans). The mass error and 1/K0 are intensity-weighted means over the
          // same points, which collapses the per-scan centroid scatter of one ion.
          double total = 0.0, weighted_ppm = 0.0, weighted_im = 0.0;
          for (auto point = spectrum.MZBegin(target - tolerance); point != spectrum.end() && point->getMZ() <= target + tolerance; ++point)
          {
            const double value = point->getIntensity();
            if (!std::isfinite(value) || value <= 0.0)
            {
              continue;
            }
            double point_im = 0.0;
            if (has_im && member.library_im > 0.0)
            {
              point_im = mobilities != nullptr ? (*mobilities)[point - spectrum.begin()] : scalar_im;
              if (!std::isfinite(point_im) || std::fabs(point_im - member.library_im) > options_.search_im)
              {
                continue;
              }
            }
            total += value;
            weighted_ppm += value * (point->getMZ() - target) / target * 1e6;
            weighted_im += value * point_im;
          }
          intensity[f] = static_cast<float>(total);
          if (total > 0.0)
          {
            ppm[f] = static_cast<float>(weighted_ppm / total);
            im[f] = static_cast<float>(weighted_im / total);
          }
        }
        if (merge)
        {
          // Different isolation windows of one frame can both cover a precursor; keep the stronger one per fragment.
          for (std::size_t f = 0; f < member.fragment_count; ++f)
          {
            if (intensity[f] > member.intensity.back()[f])
            {
              member.intensity.back()[f] = intensity[f];
              member.ppm.back()[f] = ppm[f];
              member.im.back()[f] = im[f];
            }
          }
          continue;
        }
        member.rt.push_back(rt);
        member.intensity.push_back(intensity);
        member.ppm.push_back(ppm);
        member.im.push_back(im);
      }
    }
  }

  std::vector<CalibrationCollector::Candidate> CalibrationCollector::candidates_(std::size_t index, std::size_t half_window,
                                                                                 double rt_low, double rt_high,
                                                                                 std::size_t max_candidates,
                                                                                 std::size_t min_fragments) const
  {
    const Member& member = panel_[index];
    const std::size_t n = member.rt.size();
    std::vector<double> sum(n, 0.0);
    for (std::size_t s = 0; s < n; ++s)
    {
      for (std::size_t f = 0; f < member.fragment_count; ++f)
      {
        sum[s] += member.intensity[s][f];
      }
    }
    // Local maxima of the summed trace, strongest first, at least 2 half windows apart.
    std::vector<std::size_t> maxima;
    for (std::size_t s = 1; s + 1 < n; ++s)
    {
      if (sum[s] > 0.0 && sum[s] >= sum[s - 1] && sum[s] > sum[s + 1] && member.rt[s] >= rt_low && member.rt[s] <= rt_high)
      {
        maxima.push_back(s);
      }
    }
    std::stable_sort(maxima.begin(), maxima.end(), [&](std::size_t a, std::size_t b) { return sum[a] > sum[b]; });
    std::vector<Candidate> out;
    std::vector<std::size_t> picked;
    double library_norm = 0.0;
    for (std::size_t f = 0; f < member.fragment_count; ++f)
    {
      library_norm += member.library_intensity[f];
    }
    for (const std::size_t apex : maxima)
    {
      if (out.size() >= max_candidates)
      {
        break;
      }
      if (std::any_of(picked.begin(), picked.end(), [&](std::size_t p) { return (apex > p ? apex - p : p - apex) < 2 * half_window; }))
      {
        continue;
      }
      picked.push_back(apex);
      if (apex < half_window || apex + half_window >= n)
      {
        continue;
      }
      std::vector<std::size_t> present;
      for (std::size_t f = 0; f < member.fragment_count; ++f)
      {
        if (member.intensity[apex][f] > 0.0f)
        {
          present.push_back(f);
        }
      }
      if (present.size() < min_fragments)
      {
        continue;
      }
      // Co-elution: correlation of each fragment trace with the sum of the other fragments over the
      // window; the mean of the four best fragments tolerates two interfered or missing fragments.
      const std::size_t width = 2 * half_window + 1;
      std::vector<double> total(width, 0.0);
      std::vector<std::vector<double>> traces(member.fragment_count, std::vector<double>(width, 0.0));
      for (std::size_t f = 0; f < member.fragment_count; ++f)
      {
        for (std::size_t s = 0; s < width; ++s)
        {
          traces[f][s] = member.intensity[apex - half_window + s][f];
          total[s] += traces[f][s];
        }
      }
      std::vector<double> correlations;
      for (std::size_t f = 0; f < member.fragment_count; ++f)
      {
        std::vector<double> rest(width);
        for (std::size_t s = 0; s < width; ++s)
        {
          rest[s] = total[s] - traces[f][s];
        }
        const double r = pearson_(traces[f], rest);
        if (std::isfinite(r))
        {
          correlations.push_back(r);
        }
      }
      if (correlations.size() < CORRELATED_FRAGMENTS)
      {
        continue;
      }
      std::partial_sort(correlations.begin(), correlations.begin() + CORRELATED_FRAGMENTS, correlations.end(), std::greater<>());
      const double correlation = std::accumulate(correlations.begin(), correlations.begin() + CORRELATED_FRAGMENTS, 0.0) / CORRELATED_FRAGMENTS;
      // Cosine of sqrt(apex intensities) vs sqrt(library intensities).
      double dot = 0.0, norm_x = 0.0, norm_y = 0.0;
      for (std::size_t f = 0; f < member.fragment_count; ++f)
      {
        const double x = std::sqrt(static_cast<double>(member.intensity[apex][f]));
        const double y = std::sqrt(member.library_intensity[f]);
        dot += x * y;
        norm_x += x * x;
        norm_y += y * y;
      }
      Candidate candidate;
      candidate.member = index;
      candidate.apex = apex;
      candidate.rt = member.rt[apex];
      candidate.correlation = correlation;
      candidate.cosine = norm_x > 0.0 && norm_y > 0.0 && library_norm > 0.0 ? dot / std::sqrt(norm_x * norm_y) : 0.0;
      candidate.score = candidate.correlation * candidate.cosine;
      candidate.fragments = present.size();
      out.push_back(candidate);
    }
    return out;
  }

  CalibrationResult CalibrationCollector::calibrate(CalibrationDiagnostics* diagnostics) const
  {
    if (!saw_spectrum_ || !(last_rt_ > first_rt_))
    {
      throw std::runtime_error("Calibration saw no MS2 spectra spanning a positive RT range");
    }
    if (diagnostics != nullptr)
    {
      diagnostics->rows.assign(panel_.size(), {});
      diagnostics->first_rt = first_rt_;
      diagnostics->last_rt = last_rt_;
      for (std::size_t m = 0; m < panel_.size(); ++m)
      {
        auto& row = diagnostics->rows[m];
        row.precursor_id = panel_[m].id;
        row.is_decoy = panel_[m].is_decoy;
        row.library_rt = panel_[m].library_rt;
        row.spectra = panel_[m].rt.size();
      }
    }
    // Library RT of the target a member is searched around (a decoy uses its target's).
    const auto target_library_rt = [this](std::size_t m)
    {
      return panel_[panel_[m].is_decoy ? panel_[m].partner : m].library_rt;
    };
    if (window_)
    {
      // Windowed second pass: stages 1 and 2 come from the first pass; stage 3 is centred on its stage-1 map.
      CalibrationResult result;
      for (const auto& member : panel_)
      {
        (member.is_decoy ? result.panel_decoys : result.panel_targets) += 1;
      }
      result.pairs = panel_.size() / 2;
      result.windowed = true;
      result.stage1_kept = window_->stage1_kept;
      result.stage1_sigma_s = window_->stage1_sigma_s;
      result.peak_fwhm_s = window_->peak_fwhm_s;
      result.half_window_spectra = window_->half_window_spectra;
      result.stage3_half_window_s = window_->half_window_s;
      const RTMap& map = window_->map;
      if (diagnostics != nullptr)
      {
        diagnostics->peak_fwhm_s = window_->peak_fwhm_s;
        diagnostics->half_window_spectra = window_->half_window_spectra;
        for (std::size_t m = 0; m < panel_.size(); ++m)
        {
          diagnostics->rows[m].stage1_map_rt = map(target_library_rt(m));
        }
      }
      finish_(result, [&map](double library_rt) { return map(library_rt); }, window_->half_window_spectra, diagnostics);
      return result;
    }
    CalibrationResult result;
    for (const auto& member : panel_)
    {
      (member.is_decoy ? result.panel_decoys : result.panel_targets) += 1;
    }
    result.pairs = panel_.size() / 2;
    const double run_span = last_rt_ - first_rt_;

    // Typical spectrum spacing per precursor (cycle time).
    std::vector<double> spacing;
    for (const auto& member : panel_)
    {
      for (std::size_t s = 1; s < member.rt.size(); s += std::max<std::size_t>(1, member.rt.size() / 8))
      {
        spacing.push_back(member.rt[s] - member.rt[s - 1]);
      }
    }
    const double cycle_s = std::max(1e-3, median_(spacing));

    // ---- Stage 1: whole-run anchors by paired target-decoy competition (no RT information), then the
    // per-bin RT-mode filter (DIATraceExtractor concept) as a guard against a few remaining outliers.
    std::vector<Candidate> stage1_best(panel_.size());
    std::vector<double> stage1_score(panel_.size(), -1.0);
    std::vector<bool> stage1_present(panel_.size(), false), panel_decoy(panel_.size(), false);
    std::vector<std::size_t> panel_partner(panel_.size());
    for (std::size_t m = 0; m < panel_.size(); ++m)
    {
      panel_decoy[m] = panel_[m].is_decoy;
      panel_partner[m] = panel_[m].partner;
      for (const auto& candidate : candidates_(m, options_.seed_half_window, -std::numeric_limits<double>::infinity(),
                                               std::numeric_limits<double>::infinity(), options_.candidates, options_.min_fragments))
      {
        if (!stage1_present[m] || candidate.score > stage1_score[m])
        {
          stage1_present[m] = true;
          stage1_score[m] = candidate.score;
          stage1_best[m] = candidate;
        }
      }
    }
    const auto stage1_q = pairedQValues_(stage1_score, panel_decoy, panel_partner, stage1_present);
    if (diagnostics != nullptr)
    {
      for (std::size_t m = 0; m < panel_.size(); ++m)
      {
        if (stage1_present[m])
        {
          diagnostics->rows[m].stage1_rt = stage1_best[m].rt;
          diagnostics->rows[m].stage1_score = stage1_score[m];
          diagnostics->rows[m].stage1_q = stage1_q[m];
        }
      }
    }
    std::vector<std::pair<double, double>> stage1;
    std::vector<double> stage1_fwhm;
    for (std::size_t m = 0; m < panel_.size(); ++m)
    {
      if (panel_decoy[m] || !stage1_present[m] || stage1_q[m] > options_.stage1_fdr)
      {
        continue;
      }
      const Candidate* best = &stage1_best[m];
      ++result.stage1_targets;
      stage1.emplace_back(panel_[m].library_rt, best->rt);
      // Width at half the apex of the summed trace, in spectra.
      const auto& member = panel_[m];
      const auto total = [&](std::size_t s)
      {
        return std::accumulate(member.intensity[s].begin(), member.intensity[s].end(), 0.0);
      };
      const double half = 0.5 * total(best->apex);
      std::size_t left = best->apex, right = best->apex;
      while (left > 0 && best->apex - left < 64 && total(left - 1) > half)
      {
        --left;
      }
      while (right + 1 < member.rt.size() && right - best->apex < 64 && total(right + 1) > half)
      {
        ++right;
      }
      stage1_fwhm.push_back(static_cast<double>(right - left + 1));
    }
    {
      // Decoys that won their pair above the weakest accepted stage-1 target.
      double weakest = std::numeric_limits<double>::infinity();
      for (std::size_t m = 0; m < panel_.size(); ++m)
      {
        if (!panel_decoy[m] && stage1_present[m] && stage1_q[m] <= options_.stage1_fdr)
        {
          weakest = std::min(weakest, stage1_score[m]);
        }
      }
      for (std::size_t m = 0; m < panel_.size(); ++m)
      {
        const std::size_t t = panel_partner[m];
        result.stage1_decoys += panel_decoy[m] && stage1_present[m] && (!stage1_present[t] || stage1_score[m] >= stage1_score[t]) &&
                                stage1_score[m] >= weakest;
      }
    }
    if (stage1.size() < 30)
    {
      throw std::runtime_error("Calibration found only " + std::to_string(stage1.size()) + " whole-run stage-1 anchors at q <= " +
                               std::to_string(options_.stage1_fdr) + " of " + std::to_string(result.panel_targets) +
                               " panel targets; the panel or extraction windows do not match this run");
    }
    std::sort(stage1.begin(), stage1.end());
    {
      // Per library-RT bin, the run RT with the most anchors within +- bandwidth; drop anchors far from the mode curve.
      const std::size_t bins = std::max<std::size_t>(1, std::min(options_.mode_bins, stage1.size() / 10));
      std::vector<double> mode_x, mode_y;
      for (std::size_t b = 0; b < bins; ++b)
      {
        const std::size_t lo = stage1.size() * b / bins, hi = stage1.size() * (b + 1) / bins;
        if (hi - lo < 5)
        {
          continue;
        }
        std::vector<double> y;
        for (std::size_t i = lo; i < hi; ++i)
        {
          y.push_back(stage1[i].second);
        }
        std::sort(y.begin(), y.end());
        std::size_t best = 0, best_count = 0, k = 0;
        for (std::size_t i = 0; i < y.size(); ++i)
        {
          while (y[i] - y[k] > 2.0 * options_.mode_bandwidth_s)
          {
            ++k;
          }
          if (i - k + 1 > best_count)
          {
            best_count = i - k + 1;
            best = (i + k) / 2;
          }
        }
        mode_x.push_back(stage1[(lo + hi) / 2].first);
        mode_y.push_back(y[best]);
      }
      std::vector<std::pair<double, double>> kept;
      for (const auto& point : stage1)
      {
        if (std::fabs(point.second - interpolate_(mode_x, mode_y, point.first)) <= 4.0 * options_.mode_bandwidth_s)
        {
          kept.push_back(point);
        }
      }
      stage1 = std::move(kept);
    }
    result.stage1_kept = stage1.size();
    if (stage1.size() < 20)
    {
      throw std::runtime_error("Calibration stage-1 anchors do not agree on an RT curve");
    }
    const RTMap map1 = fitRTMap(stage1, mapBins_(stage1.size()));
    result.stage1_sigma_s = map1.sigma;
    if (diagnostics != nullptr)
    {
      for (std::size_t m = 0; m < panel_.size(); ++m)
      {
        diagnostics->rows[m].stage1_map_rt = map1(target_library_rt(m));
      }
    }
    const std::size_t half_window = std::clamp<std::size_t>(static_cast<std::size_t>(std::lround(median_(stage1_fwhm))), 3, 64);
    result.half_window_spectra = half_window;
    result.peak_fwhm_s = median_(stage1_fwhm) * cycle_s;
    if (diagnostics != nullptr)
    {
      diagnostics->cycle_s = cycle_s;
      diagnostics->peak_fwhm_s = result.peak_fwhm_s;
      diagnostics->half_window_spectra = half_window;
    }

    // ---- Stage 2: identification-free latent peak/no-peak mixture over all target candidates.
    // Library RT normalized over the panel's target span; run RT over the measured run interval.
    double lib_lo = std::numeric_limits<double>::infinity(), lib_hi = -lib_lo;
    for (const auto& member : panel_)
    {
      if (!member.is_decoy)
      {
        lib_lo = std::min(lib_lo, member.library_rt);
        lib_hi = std::max(lib_hi, member.library_rt);
      }
    }
    if (!(lib_hi > lib_lo))
    {
      throw std::runtime_error("Calibration panel has no library RT spread");
    }
    struct Hypothesis
    {
      double x = 0.0;
      std::vector<double> y, contrast;
    };
    std::vector<Hypothesis> hypotheses;
    std::vector<double> deviations;
    std::vector<std::vector<double>> raw_scores;
    for (std::size_t m = 0; m < panel_.size(); ++m)
    {
      if (panel_[m].is_decoy)
      {
        continue;
      }
      const auto candidates = candidates_(m, half_window, -std::numeric_limits<double>::infinity(),
                                          std::numeric_limits<double>::infinity(), options_.candidates, options_.min_fragments);
      if (candidates.empty())
      {
        continue;
      }
      Hypothesis h;
      h.x = (panel_[m].library_rt - lib_lo) / (lib_hi - lib_lo);
      std::vector<double> scores;
      for (const auto& candidate : candidates)
      {
        h.y.push_back((candidate.rt - first_rt_) / run_span);
        scores.push_back(candidate.score);
      }
      const double centre = median_(scores);
      for (double& score : scores)
      {
        score -= centre;
        deviations.push_back(std::fabs(score));
      }
      raw_scores.push_back(std::move(scores));
      hypotheses.push_back(std::move(h));
    }
    if (hypotheses.size() < 30)
    {
      throw std::runtime_error("Calibration stage 2 has too few candidate-bearing targets");
    }
    // One panel-wide robust scale for the contrast, floored at 0.05; contrast clipped to [-3, 5].
    const double scale = std::max(0.05, 1.4826 * median_(deviations));
    for (std::size_t i = 0; i < hypotheses.size(); ++i)
    {
      for (double score : raw_scores[i])
      {
        hypotheses[i].contrast.push_back(std::clamp(score / scale, -3.0, 5.0));
      }
    }
    // Knots at library-RT quantiles of the hypotheses; initialized from the stage-1 map.
    // v0.13: equally spaced knots on normalized library RT; the prior heights come from the stage-1 map.
    std::vector<double> knot_x, knot_y;
    for (std::size_t k = 0; k < options_.em_knots; ++k)
    {
      knot_x.push_back(static_cast<double>(k) / static_cast<double>(options_.em_knots - 1));
    }
    for (double x : knot_x)
    {
      knot_y.push_back((map1(lib_lo + x * (lib_hi - lib_lo)) - first_rt_) / run_span);
    }
    // v0.13 keeps sigma and the presence prior fixed; they are not re-estimated.
    const double sigma = options_.em_sigma_fraction;
    const double presence = options_.presence_prior;
    const std::vector<double> prior_y = knot_y;
    const std::size_t knots = knot_x.size();
    const double hypothesis_count = static_cast<double>(hypotheses.size());
    for (result.em_iterations = 0; result.em_iterations < options_.em_max_iterations; ++result.em_iterations)
    {
      // E step: responsibilities of each candidate peak and of the explicit no-peak hypothesis
      // (uniform density 1 on normalized run time), with the K_i correction.
      // M step: minimize mean over hypotheses of sum_j r_ij (y_ij - f(x_i))^2 / (2 sigma^2)
      // + prior_penalty * mean_k (theta_k - prior_k)^2 + curve_penalty * mean_k (second difference)^2,
      // a linear system in the knot heights theta (hat basis), then enforce a minimum knot step.
      std::vector<double> a(knots * knots, 0.0), b(knots, 0.0);
      double responsibility_sum = 0.0;
      for (const auto& h : hypotheses)
      {
        const double centre = interpolate_(knot_x, knot_y, h.x);
        const double k = static_cast<double>(h.y.size());
        std::vector<double> r(h.y.size());
        double total = 1.0 - presence;
        for (std::size_t j = 0; j < h.y.size(); ++j)
        {
          const double z = (h.y[j] - centre) / sigma;
          r[j] = presence / k * std::exp(options_.score_beta * h.contrast[j]) * std::exp(-0.5 * z * z) / (sigma * std::sqrt(2.0 * PI));
          total += r[j];
        }
        const std::size_t upper = std::min<std::size_t>(knots - 1,
          std::max<std::size_t>(1, std::upper_bound(knot_x.begin(), knot_x.end(), h.x) - knot_x.begin()));
        const double t = std::clamp((h.x - knot_x[upper - 1]) / (knot_x[upper] - knot_x[upper - 1]), 0.0, 1.0);
        const std::size_t lo = upper - 1;
        const double w_lo = 1.0 - t, w_hi = t;
        for (std::size_t j = 0; j < h.y.size(); ++j)
        {
          const double weight = r[j] / total / (sigma * sigma * hypothesis_count);
          a[lo * knots + lo] += weight * w_lo * w_lo;
          a[lo * knots + upper] += weight * w_lo * w_hi;
          a[upper * knots + lo] += weight * w_lo * w_hi;
          a[upper * knots + upper] += weight * w_hi * w_hi;
          b[lo] += weight * w_lo * h.y[j];
          b[upper] += weight * w_hi * h.y[j];
          responsibility_sum += r[j] / total;
        }
      }
      // Normal equations: the data term's gradient carries 1/sigma^2, each penalty lambda * mean(.)^2 carries 2 lambda / count.
      const double prior_penalty = 2.0 * 0.12 / static_cast<double>(knots);
      const double curve_penalty = 2.0 * 0.2 / static_cast<double>(std::max<std::size_t>(1, knots - 2));
      for (std::size_t k = 0; k < knots; ++k)
      {
        a[k * knots + k] += prior_penalty;
        b[k] += prior_penalty * prior_y[k];
      }
      for (std::size_t k = 1; k + 1 < knots; ++k)
      {
        const std::array<std::size_t, 3> idx{k - 1, k, k + 1};
        const std::array<double, 3> coef{1.0, -2.0, 1.0};
        for (std::size_t u = 0; u < 3; ++u)
        {
          for (std::size_t v = 0; v < 3; ++v)
          {
            a[idx[u] * knots + idx[v]] += curve_penalty * coef[u] * coef[v];
          }
        }
      }
      // Gaussian elimination with partial pivoting on the small dense system.
      std::vector<double> next = b;
      for (std::size_t col = 0; col < knots; ++col)
      {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < knots; ++row)
        {
          if (std::fabs(a[row * knots + col]) > std::fabs(a[pivot * knots + col]))
          {
            pivot = row;
          }
        }
        for (std::size_t c = 0; c < knots; ++c)
        {
          std::swap(a[col * knots + c], a[pivot * knots + c]);
        }
        std::swap(next[col], next[pivot]);
        for (std::size_t row = col + 1; row < knots; ++row)
        {
          const double factor = a[row * knots + col] / a[col * knots + col];
          for (std::size_t c = col; c < knots; ++c)
          {
            a[row * knots + c] -= factor * a[col * knots + c];
          }
          next[row] -= factor * next[col];
        }
      }
      for (std::size_t col = knots; col-- > 0;)
      {
        for (std::size_t c = col + 1; c < knots; ++c)
        {
          next[col] -= a[col * knots + c] * next[c];
        }
        next[col] /= a[col * knots + col];
      }
      // Monotone with a minimum normalized step of 0.002: isotonic fit of theta_k - 0.002 k.
      for (std::size_t k = 0; k < knots; ++k)
      {
        next[k] -= 0.002 * static_cast<double>(k);
      }
      next = isotonic_(next, std::vector<double>(knots, 1.0));
      for (std::size_t k = 0; k < knots; ++k)
      {
        next[k] = std::clamp(next[k] + 0.002 * static_cast<double>(k), 0.0, 1.0);
      }
      double change = 0.0;
      for (std::size_t k = 0; k < knots; ++k)
      {
        change = std::max(change, std::fabs(next[k] - knot_y[k]));
      }
      knot_y = std::move(next);
      result.em_presence = responsibility_sum / hypothesis_count; // Mean peak responsibility (diagnostic).
      if (change < 1e-5)
      {
        result.em_converged = true;
        ++result.em_iterations;
        break;
      }
    }
    result.em_sigma_s = sigma * run_span;
    for (std::size_t k = 0; k < knot_x.size(); ++k)
    {
      result.em_knots_library.push_back(lib_lo + knot_x[k] * (lib_hi - lib_lo));
      result.em_knots_run.push_back(first_rt_ + knot_y[k] * run_span);
    }
    const auto em_map = [&](double library_rt)
    {
      return interpolate_(result.em_knots_library, result.em_knots_run, library_rt);
    };
    if (diagnostics != nullptr)
    {
      for (std::size_t m = 0; m < panel_.size(); ++m)
      {
        diagnostics->rows[m].em_rt = em_map(target_library_rt(m));
      }
    }

    // ---- Stage 3: paired target-decoy anchors inside the stage-2 window.
    // The EM sigma is a fixed model width, not a measured spread; the window uses the stage-1 anchor spread.
    result.stage3_half_window_s = std::max(options_.rt_sigma_multiple * result.stage1_sigma_s, options_.min_window_s);
    try
    {
      finish_(result, em_map, half_window, diagnostics);
    }
    catch (const std::runtime_error&)
    {
      if (diagnostics != nullptr)
      {
        // Stage 1 worked; keep its map so a caller may try a windowed second pass.
        CalibrationWindow window;
        window.map = map1;
        window.half_window_s = result.stage3_half_window_s;
        window.half_window_spectra = half_window;
        window.stage1_kept = result.stage1_kept;
        window.stage1_sigma_s = result.stage1_sigma_s;
        window.peak_fwhm_s = result.peak_fwhm_s;
        diagnostics->window = window;
      }
      throw;
    }
    return result;
  }

  void CalibrationCollector::finish_(CalibrationResult& result, const std::function<double(double)>& centre_of,
                                     std::size_t half_window, CalibrationDiagnostics* diagnostics) const
  {
    const auto target_library_rt = [this](std::size_t m)
    {
      return panel_[panel_[m].is_decoy ? panel_[m].partner : m].library_rt;
    };
    std::vector<double> score(panel_.size(), -1.0);
    std::vector<bool> present(panel_.size(), false), is_decoy(panel_.size(), false);
    std::vector<std::size_t> partner(panel_.size());
    std::vector<Candidate> best(panel_.size());
    for (std::size_t m = 0; m < panel_.size(); ++m)
    {
      is_decoy[m] = panel_[m].is_decoy;
      partner[m] = panel_[m].partner;
      // A decoy is searched around its target's predicted RT, like the target.
      const double centre = centre_of(target_library_rt(m));
      if (diagnostics != nullptr)
      {
        diagnostics->rows[m].stage3_centre_rt = centre;
      }
      const auto candidates = candidates_(m, half_window, centre - result.stage3_half_window_s, centre + result.stage3_half_window_s,
                                          options_.candidates, options_.min_fragments);
      for (const auto& candidate : candidates)
      {
        if (!present[m] || candidate.score > score[m])
        {
          present[m] = true;
          score[m] = candidate.score;
          best[m] = candidate;
        }
      }
    }
    const auto q = pairedQValues_(score, is_decoy, partner, present);
    if (diagnostics != nullptr)
    {
      for (std::size_t m = 0; m < panel_.size(); ++m)
      {
        if (present[m])
        {
          diagnostics->rows[m].stage3_rt = best[m].rt;
          diagnostics->rows[m].stage3_score = score[m];
          diagnostics->rows[m].stage3_q = q[m];
        }
      }
    }
    std::vector<std::pair<double, double>> anchor_points;
    for (std::size_t m = 0; m < panel_.size(); ++m)
    {
      if (!present[m])
      {
        continue;
      }
      if (is_decoy[m])
      {
        continue;
      }
      if (q[m] > options_.anchor_fdr)
      {
        continue;
      }
      const auto& member = panel_[m];
      CalibrationAnchor anchor;
      anchor.precursor_id = member.id;
      anchor.library_rt = member.library_rt;
      anchor.run_rt = best[m].rt;
      anchor.library_im = member.library_im;
      anchor.correlation = best[m].correlation;
      anchor.cosine = best[m].cosine;
      anchor.q_value = q[m];
      anchor.fragments = best[m].fragments;
      std::vector<double> ims;
      for (std::size_t f = 0; f < member.fragment_count; ++f)
      {
        if (member.intensity[best[m].apex][f] > 0.0f && member.im[best[m].apex][f] > 0.0f)
        {
          ims.push_back(member.im[best[m].apex][f]);
        }
      }
      anchor.observed_im = median_(ims);
      result.anchors.push_back(anchor);
      anchor_points.emplace_back(anchor.library_rt, anchor.run_rt);
    }
    result.anchor_targets = result.anchors.size();
    // Decoys that won their pair at or above the weakest accepted anchor score.
    double weakest = std::numeric_limits<double>::infinity();
    for (std::size_t m = 0; m < panel_.size(); ++m)
    {
      if (!is_decoy[m] && present[m] && q[m] <= options_.anchor_fdr)
      {
        weakest = std::min(weakest, score[m]);
      }
    }
    for (std::size_t m = 0; m < panel_.size(); ++m)
    {
      const std::size_t t = partner[m];
      result.anchor_decoys += is_decoy[m] && present[m] && (!present[t] || score[m] >= score[t]) && score[m] >= weakest;
    }
    if (anchor_points.size() < 20)
    {
      std::vector<std::pair<double, bool>> ranked;
      for (std::size_t m = 0; m < panel_.size(); ++m)
      {
        if (present[m])
        {
          ranked.emplace_back(score[m], is_decoy[m]);
        }
      }
      std::sort(ranked.rbegin(), ranked.rend());
      std::size_t top_decoys = 0;
      for (std::size_t i = 0; i < std::min<std::size_t>(100, ranked.size()); ++i)
      {
        top_decoys += ranked[i].second;
      }
      throw std::runtime_error("Calibration found only " + std::to_string(anchor_points.size()) + " anchors at the anchor FDR (stage1 kept " +
                               std::to_string(result.stage1_kept) + ", stage1 sigma " + std::to_string(result.stage1_sigma_s) +
                               " s, EM peak responsibility " + std::to_string(result.em_presence) +
                               ", stage3 half window " + std::to_string(result.stage3_half_window_s) + " s, " +
                               (result.windowed ? std::string("windowed pass, ") : std::string()) + "scored " +
                               std::to_string(ranked.size()) + ", decoys in top 100 " + std::to_string(top_decoys) + ")");
    }

    // ---- Stage 4: final RT map and tolerances, fragment mass and mobility models on the anchors.
    result.rt_map = fitRTMap(anchor_points, mapBins_(anchor_points.size()));
    const double cut = std::max(3.0 * result.rt_map.sigma, 3.0);
    std::vector<double> rt_residuals;
    std::unordered_map<std::string, std::size_t> member_of;
    for (std::size_t m = 0; m < panel_.size(); ++m)
    {
      member_of.emplace(panel_[m].id, m);
    }
    std::vector<double> ppm_values;
    std::vector<std::pair<double, double>> im_pairs;
    for (auto& anchor : result.anchors)
    {
      const double residual = anchor.run_rt - result.rt_map(anchor.library_rt);
      anchor.used = std::fabs(residual) <= cut;
      if (!anchor.used)
      {
        continue;
      }
      rt_residuals.push_back(std::fabs(residual));
      const auto& member = panel_[member_of.at(anchor.precursor_id)];
      const auto& candidate = best[member_of.at(anchor.precursor_id)];
      for (std::size_t f = 0; f < member.fragment_count; ++f)
      {
        if (member.intensity[candidate.apex][f] > 0.0f)
        {
          ppm_values.push_back(member.ppm[candidate.apex][f]);
        }
      }
      if (anchor.observed_im > 0.0 && anchor.library_im > 0.0)
      {
        im_pairs.emplace_back(anchor.library_im, anchor.observed_im);
      }
    }
    result.rt_p50_s = quantile_(rt_residuals, 0.50);
    result.rt_p95_s = quantile_(rt_residuals, 0.95);
    result.rt_p99_s = quantile_(rt_residuals, 0.99);
    result.rt_tolerance_s = options_.rt_sigma_multiple * result.rt_map.sigma;
    result.ms2_points = ppm_values.size();
    if (ppm_values.size() >= 50)
    {
      result.ms2_offset_ppm = median_(ppm_values);
      result.ms2_sigma_ppm = robustSigma_(ppm_values, result.ms2_offset_ppm);
      result.ms2_tolerance_ppm = std::max(options_.ms2_sigma_multiple * result.ms2_sigma_ppm, 0.1);
    }
    result.im_points = im_pairs.size();
    if (saw_im_ && im_pairs.size() >= 30)
    {
      // Robust line: least squares on the core (|residual| < 4 sigma about the median offset), twice.
      std::vector<double> offsets;
      for (const auto& [library_im, observed_im] : im_pairs)
      {
        offsets.push_back(observed_im - library_im);
      }
      double intercept = median_(offsets), slope = 1.0;
      double im_sigma = std::max(robustSigma_(offsets, intercept), 1e-4);
      for (int pass = 0; pass < 2; ++pass)
      {
        double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0, n = 0.0;
        for (const auto& [x, y] : im_pairs)
        {
          if (std::fabs(y - (intercept + slope * x)) >= 4.0 * im_sigma)
          {
            continue;
          }
          sx += x;
          sy += y;
          sxx += x * x;
          sxy += x * y;
          n += 1.0;
        }
        const double det = n * sxx - sx * sx;
        if (n < 30 || det <= 0.0)
        {
          break;
        }
        slope = (n * sxy - sx * sy) / det;
        intercept = (sy - slope * sx) / n;
        std::vector<double> residuals;
        for (const auto& [x, y] : im_pairs)
        {
          residuals.push_back(y - (intercept + slope * x));
        }
        im_sigma = std::max(robustSigma_(residuals, 0.0), 1e-4);
      }
      result.has_im = true;
      result.im_intercept = intercept;
      result.im_slope = slope;
      result.im_sigma = im_sigma;
      result.im_tolerance = options_.im_sigma_multiple * im_sigma;
    }
  }

  void writeCalibrationJSON(const CalibrationResult& r, const CalibrationOptions& o, std::ostream& out)
  {
    const auto list = [&](const std::vector<double>& values)
    {
      out << '[';
      for (std::size_t i = 0; i < values.size(); ++i)
      {
        out << (i ? "," : "") << values[i];
      }
      out << ']';
    };
    out << std::setprecision(10) << "{\n"
        << "  \"version\": 1,\n  \"experimental\": true,\n"
        << "  \"method\": \"co-elution anchors + identification-free latent peak/no-peak RT mixture + paired target-decoy anchors\",\n"
        << "  \"options\": {\"panel_targets\": " << o.panel_targets << ", \"search_ppm\": " << o.search_ppm
        << ", \"search_im\": " << o.search_im << ", \"candidates\": " << o.candidates
        << ", \"stage1_fdr\": " << o.stage1_fdr << ", \"correlated_fragments\": " << CORRELATED_FRAGMENTS
        << ", \"presence_prior\": " << o.presence_prior << ", \"em_sigma_fraction\": " << o.em_sigma_fraction << ", \"score_beta\": " << o.score_beta << ", \"em_knots\": " << o.em_knots
        << ", \"anchor_fdr\": " << o.anchor_fdr << ", \"rt_sigma_multiple\": " << o.rt_sigma_multiple
        << ", \"ms2_sigma_multiple\": " << o.ms2_sigma_multiple << ", \"im_sigma_multiple\": " << o.im_sigma_multiple << "},\n"
        << "  \"panel\": {\"targets\": " << r.panel_targets << ", \"decoys\": " << r.panel_decoys << "},\n"
        << "  \"stage1\": {\"targets\": " << r.stage1_targets << ", \"decoys_above_cut\": " << r.stage1_decoys
        << ", \"mode_kept\": " << r.stage1_kept << ", \"sigma_s\": " << r.stage1_sigma_s
        << ", \"peak_fwhm_s\": " << r.peak_fwhm_s << ", \"half_window_spectra\": " << r.half_window_spectra << "},\n"
        << "  \"stage2\": {\"iterations\": " << r.em_iterations << ", \"converged\": " << (r.em_converged ? "true" : "false")
        << ", \"sigma_s\": " << r.em_sigma_s << ", \"mean_peak_responsibility\": " << r.em_presence << ", \"knots_library_rt\": ";
    list(r.em_knots_library);
    out << ", \"knots_run_rt_s\": ";
    list(r.em_knots_run);
    out << "},\n  \"stage3\": {\"half_window_s\": " << r.stage3_half_window_s << ", \"anchor_targets\": " << r.anchor_targets
        << ", \"decoys_above_cut\": " << r.anchor_decoys;
    if (r.windowed)
    {
      // Only written for the experimental windowed pass, so the default file is unchanged.
      out << ", \"windowed_second_pass\": true, \"centre\": \"stage1_map\"";
    }
    out << "},\n"
        << "  \"rt\": {\"knots_library_rt\": ";
    list(r.rt_map.library_rt);
    out << ", \"knots_run_rt_s\": ";
    list(r.rt_map.run_rt);
    out << ", \"anchors_used\": " << r.rt_map.anchors << ", \"sigma_s\": " << r.rt_map.sigma << ", \"abs_residual_p50_s\": " << r.rt_p50_s
        << ", \"abs_residual_p95_s\": " << r.rt_p95_s << ", \"abs_residual_p99_s\": " << r.rt_p99_s
        << ", \"tolerance_half_width_s\": " << r.rt_tolerance_s << "},\n"
        << "  \"ms2_mz\": {\"offset_ppm\": " << r.ms2_offset_ppm << ", \"sigma_ppm\": " << r.ms2_sigma_ppm
        << ", \"tolerance_half_width_ppm\": " << r.ms2_tolerance_ppm << ", \"points\": " << r.ms2_points << "},\n"
        << "  \"im\": {\"available\": " << (r.has_im ? "true" : "false") << ", \"intercept_vssc\": " << r.im_intercept
        << ", \"slope\": " << r.im_slope << ", \"sigma_vssc\": " << r.im_sigma << ", \"tolerance_half_width_vssc\": " << r.im_tolerance
        << ", \"points\": " << r.im_points << "}\n}\n";
  }

  void writeCalibrationDiagnostics(const CalibrationDiagnostics& d, std::ostream& out)
  {
    out << "# first_rt_s=" << std::setprecision(10) << d.first_rt << " last_rt_s=" << d.last_rt << " cycle_s=" << d.cycle_s
        << " peak_fwhm_s=" << d.peak_fwhm_s << " half_window_spectra=" << d.half_window_spectra
        << " failure=" << (d.failure.empty() ? std::string("none") : d.failure) << '\n'
        << "precursor_id\tdecoy\tlibrary_rt\tspectra\tstage1_rt\tstage1_score\tstage1_q\tstage1_map_rt\tem_rt\tstage3_centre_rt\tstage3_rt\tstage3_score\tstage3_q\n";
    for (const auto& r : d.rows)
    {
      out << r.precursor_id << '\t' << r.is_decoy << '\t' << r.library_rt << '\t' << r.spectra << '\t' << r.stage1_rt << '\t'
          << r.stage1_score << '\t' << r.stage1_q << '\t' << r.stage1_map_rt << '\t' << r.em_rt << '\t' << r.stage3_centre_rt << '\t'
          << r.stage3_rt << '\t' << r.stage3_score << '\t' << r.stage3_q << '\n';
    }
  }

  void writeCalibrationAnchors(const CalibrationResult& r, std::ostream& out)
  {
    out << "precursor_id\tlibrary_rt\trun_rt\tlibrary_im\tobserved_im\tcorrelation\tcosine\tfragments\tq_value\tused\n"
        << std::setprecision(10);
    for (const auto& a : r.anchors)
    {
      out << a.precursor_id << '\t' << a.library_rt << '\t' << a.run_rt << '\t' << a.library_im << '\t' << a.observed_im << '\t'
          << a.correlation << '\t' << a.cosine << '\t' << a.fragments << '\t' << a.q_value << '\t' << a.used << '\n';
    }
  }

  CalibrationRT readCalibrationRT(const std::string& path)
  {
    std::ifstream in(path);
    if (!in)
    {
      throw std::invalid_argument("Cannot read calibration file " + path);
    }
    const nlohmann::json json = nlohmann::json::parse(in);
    if (json.at("version").get<int>() != 1)
    {
      throw std::invalid_argument("Unsupported calibration file version in " + path);
    }
    const auto& rt = json.at("rt");
    CalibrationRT result;
    result.map.library_rt = rt.at("knots_library_rt").get<std::vector<double>>();
    result.map.run_rt = rt.at("knots_run_rt_s").get<std::vector<double>>();
    result.map.sigma = rt.at("sigma_s").get<double>();
    result.map.anchors = rt.at("anchors_used").get<std::size_t>();
    result.tolerance_s = rt.at("tolerance_half_width_s").get<double>();
    const auto& x = result.map.library_rt;
    const auto& y = result.map.run_rt;
    bool valid = x.size() >= 2 && x.size() == y.size() && std::isfinite(result.tolerance_s) && result.tolerance_s > 0.0;
    for (std::size_t i = 0; valid && i < x.size(); ++i)
    {
      valid = std::isfinite(x[i]) && std::isfinite(y[i]) && (i == 0 || x[i] >= x[i - 1]);
    }
    if (!valid)
    {
      throw std::invalid_argument("Calibration file " + path + " has an invalid RT map or tolerance");
    }
    const auto& ms2 = json.at("ms2_mz");
    result.ms2_offset_ppm = ms2.at("offset_ppm").get<double>();
    result.ms2_tolerance_ppm = ms2.at("tolerance_half_width_ppm").get<double>();
    result.has_ms2 = std::isfinite(result.ms2_tolerance_ppm) && result.ms2_tolerance_ppm > 0.0 && std::isfinite(result.ms2_offset_ppm);
    const auto& im = json.at("im");
    result.has_im = im.at("available").get<bool>();
    result.im_intercept = im.at("intercept_vssc").get<double>();
    result.im_slope = im.at("slope").get<double>();
    result.im_tolerance = im.at("tolerance_half_width_vssc").get<double>();
    if (result.has_im && !(std::isfinite(result.im_intercept) && std::isfinite(result.im_slope) && result.im_slope > 0.0 &&
                           std::isfinite(result.im_tolerance) && result.im_tolerance > 0.0))
    {
      throw std::invalid_argument("Calibration file " + path + " has an invalid 1/K0 model");
    }
    return result;
  }
}
