// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIATraces.h"

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <set>

namespace OpenMS::OpenDIA
{
  namespace
  {
    void require(bool valid, const std::string& message)
    {
      if (!valid)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
      }
    }

    bool safeText(const std::string& text)
    {
      return !text.empty() && text.find_first_of("\t\r\n") == std::string::npos;
    }
  }

  std::vector<TraceWindow> readTraceWindows(const std::string& path, std::size_t max_samples)
  {
    std::ifstream input(path);
    require(input.good(), "Could not open candidate windows: " + path);
    std::string line;
    require(static_cast<bool>(std::getline(input, line)), "Candidate window table is empty.");
    if (!line.empty() && line.back() == '\r')
    {
      line.pop_back();
    }
    std::vector<std::string> columns;
    StringUtils::split(line, '\t', columns);
    const std::array<std::string, 5> names = {"sample_id", "precursor_id", "group_id", "label", "rt_apex"};
    std::array<std::size_t, 5> indices{};
    for (std::size_t i = 0; i < names.size(); ++i)
    {
      require(std::count(columns.begin(), columns.end(), names[i]) == 1, "Expected one column: " + names[i]);
      indices[i] = std::find(columns.begin(), columns.end(), names[i]) - columns.begin();
    }
    std::vector<TraceWindow> windows;
    while (std::getline(input, line))
    {
      if (!line.empty() && line.back() == '\r')
      {
        line.pop_back();
      }
      require(!line.empty(), "Blank candidate window row.");
      std::vector<std::string> values;
      StringUtils::split(line, '\t', values);
      require(values.size() == columns.size(), "Candidate window column count differs from header.");
      require(windows.size() < max_samples, "Candidate windows exceed max_samples. Select a representative subset explicitly.");
      windows.push_back({values[indices[0]], values[indices[1]], values[indices[2]],
                         StringUtils::toInt32(values[indices[3]]), StringUtils::toDouble(values[indices[4]])});
    }
    require(input.eof(), "Error reading candidate windows.");
    return windows;
  }

  TraceCollector::TraceCollector(std::vector<TraceWindow> windows, const OpenDIA::Library& library,
                                 double mz_tolerance, double rt_half_width, double im_half_width,
                                 std::size_t max_samples, std::size_t max_tensor_bytes, CloseCallback on_close) :
    mz_tolerance_(mz_tolerance), rt_half_width_(rt_half_width), im_half_width_(im_half_width), on_close_(std::move(on_close))
  {
    require(std::isfinite(mz_tolerance_) && mz_tolerance_ > 0.0 &&
            std::isfinite(rt_half_width_) && rt_half_width_ > 0.0 &&
            std::isfinite(im_half_width_) && im_half_width_ > 0.0, "Trace tolerances must be finite and positive.");
    // ponytail: sampled tensors are retained under an explicit byte budget.
    // Disk-backed tiled accumulation is only needed if a run cannot fit in RAM.
    require(!windows.empty() && windows.size() <= max_samples, "Trace sample count is empty or exceeds max_samples.");
    std::map<std::string, std::vector<const OpenDIA::LibraryTransition*>> transitions;
    std::map<std::string, double> mobilities;
    std::map<std::string, int> charges;
    for (const auto& compound : library.getCompounds())
    {
      mobilities[compound.id] = compound.drift_time;
      charges[compound.id] = compound.charge;
    }
    for (const auto& transition : library.getTransitions())
    {
      transitions[transition.peptide_ref].push_back(&transition);
    }
    std::set<std::string> ids;
    for (auto& window : windows)
    {
      require(safeText(window.sample_id) && safeText(window.precursor_id) && safeText(window.group_id) &&
              ids.insert(window.sample_id).second && std::isfinite(window.apex_rt) && window.apex_rt >= 0.0 &&
              window.label >= -1 && window.label <= 1, "Invalid/duplicate candidate sample, group, label or RT.");
      const auto found = transitions.find(window.precursor_id);
      require(found != transitions.end() && !found->second.empty(), "Unknown precursor: " + window.precursor_id);
      auto& fragments = found->second;
      for (const auto* fragment : fragments)
      {
        require(std::isfinite(fragment->product_mz) && fragment->product_mz > 0.0 &&
                std::isfinite(fragment->precursor_mz) && fragment->precursor_mz > 0.0 &&
                std::isfinite(fragment->library_intensity) && fragment->library_intensity >= 0.0 &&
                fragment->precursor_mz == fragments.front()->precursor_mz &&
                fragment->getDecoy() == fragments.front()->getDecoy(), "Invalid or inconsistent library fragment.");
      }
      std::stable_sort(fragments.begin(), fragments.end(), [](const auto* a, const auto* b)
      {
        if (a->library_intensity != b->library_intensity)
        {
          return a->library_intensity > b->library_intensity;
        }
        return a->product_mz < b->product_mz;
      });
      Sample sample;
      sample.window = std::move(window);
      sample.precursor_mz = fragments.front()->precursor_mz;
      sample.is_decoy = fragments.front()->getDecoy();
      require(!sample.is_decoy || sample.window.label != 1, "A decoy cannot be a positive seed.");
      const auto mobility = mobilities.find(sample.window.precursor_id);
      if (mobility != mobilities.end() && std::isfinite(mobility->second) && mobility->second > 0.0)
      {
        sample.precursor_im = mobility->second;
      }
      const auto charge = charges.find(sample.window.precursor_id);
      sample.precursor_charge = charge != charges.end() ? charge->second : 0;
      for (std::size_t i = 0; i < std::min(FRAGMENTS, fragments.size()); ++i)
      {
        sample.fragment_mz[i] = fragments[i]->product_mz;
        sample.library_intensity[i] = fragments[i]->library_intensity;
        sample.fragment_charge[i] = fragments[i]->fragment_charge;
      }
      samples_.push_back(std::move(sample));
    }
    // Without any library mobility every point lands in the central mobility bin (compatible_im is never true in
    // consumeSpectrum), so storing only that bin gives the same tensor at 1/8 of the memory.
    im_bins_ = std::any_of(samples_.begin(), samples_.end(), [](const Sample& s) { return s.precursor_im > 0.0; }) ? IM_BINS : 1;
    const std::size_t stored = SAMPLE_SIZE / IM_BINS * im_bins_;
    stride_ = stored;
    if (on_close_)
    {
      // Streaming: tensors live in reusable slots between a window's first spectrum and its closing.
      max_held_ = max_tensor_bytes / (stored * sizeof(float));
      require(max_held_ > 0, "Trace tensor budget holds no window.");
      slot_of_.assign(samples_.size(), NO_SLOT);
      closed_.assign(samples_.size(), 0);
    }
    else
    {
      require(samples_.size() <= max_tensor_bytes / (stored * sizeof(float)), "Trace sample count exceeds the raw tensor byte budget.");
      intensities_.resize(samples_.size() * stored, 0.0f);
    }
    by_rt_.resize(samples_.size());
    std::iota(by_rt_.begin(), by_rt_.end(), 0);
    std::stable_sort(by_rt_.begin(), by_rt_.end(), [this](std::size_t a, std::size_t b)
    { return samples_[a].window.apex_rt < samples_[b].window.apex_rt; });
    // MS2 isolation lookup: a spectrum only visits the buckets its isolation windows overlap (see consumeSpectrum).
    double low = std::numeric_limits<double>::infinity(), high = -low;
    for (const auto& s : samples_) { low = std::min(low, s.precursor_mz); high = std::max(high, s.precursor_mz); }
    require(std::isfinite(low) && std::isfinite(high), "Trace precursor m/z must be finite.");
    require(samples_.size() <= std::numeric_limits<std::uint32_t>::max(), "Trace sample count exceeds 32-bit indices.");
    mz_bucket_origin_ = std::floor(low);
    by_mz_rt_.assign(static_cast<std::size_t>(std::floor((high - mz_bucket_origin_) / MZ_BUCKET_TH)) + 1, {});
    for (const std::size_t i : by_rt_)
    {
      by_mz_rt_[static_cast<std::size_t>(std::floor((samples_[i].precursor_mz - mz_bucket_origin_) / MZ_BUCKET_TH))].push_back(static_cast<std::uint32_t>(i));
    }
  }

  void TraceCollector::enableEvidence(double ms1_ppm, double core_half_width_s)
  {
    require(std::isfinite(ms1_ppm) && ms1_ppm > 0.0 && std::isfinite(core_half_width_s) && core_half_width_s > 0.0,
            "Evidence tolerances must be finite and positive.");
    ms1_ppm_ = ms1_ppm;
    core_half_width_ = core_half_width_s;
    evidence_.assign(samples_.size(), TraceEvidence{});
  }

  void TraceCollector::enablePeakGroupEvidence(double ppm)
  {
    require(!evidence_.empty(), "Peak-group evidence requires enableEvidence() first.");
    require(std::isfinite(ppm) && ppm > 0.0, "Peak-group evidence tolerance must be finite and positive.");
    pg_ppm_ = ppm;
    pg_evidence_.assign(samples_.size(), PeakGroupEvidence{});
  }

  std::size_t TraceCollector::rtBin_(double rt, double apex_rt) const
  {
    return std::min(RT_BINS - 1, static_cast<std::size_t>(
      std::max(0.0, (rt - apex_rt + rt_half_width_) / (2.0 * rt_half_width_)) * RT_BINS));
  }

  void TraceCollector::enableByIonCounts(const OpenDIA::Library& library, double ppm, double intensity_min)
  {
    require(!evidence_.empty(), "b/y ion counts require enableEvidence() first.");
    require(std::isfinite(ppm) && ppm > 0.0 && std::isfinite(intensity_min) && intensity_min >= 0.0,
            "b/y ion tolerance must be positive and the intensity floor non-negative.");
    byion_ppm_ = ppm;
    byion_intensity_min_ = intensity_min;
    std::map<std::string, const std::string*> sequences;
    for (const auto& compound : library.getCompounds())
    {
      sequences[compound.id] = &compound.sequence;
    }
    // Ions are generated once per precursor; windows of the same precursor share the range.
    std::map<std::string, std::array<std::uint32_t, 3>> ranges;
    byion_mz_.clear();
    byion_range_.assign(samples_.size(), {0, 0, 0});
    for (std::size_t i = 0; i < samples_.size(); ++i)
    {
      const std::string& id = samples_[i].window.precursor_id;
      auto found = ranges.find(id);
      if (found == ranges.end())
      {
        const auto sequence = sequences.find(id);
        require(sequence != sequences.end() && !sequence->second->empty(), "No library sequence for precursor: " + id);
        AASequence peptide;
        try
        {
          peptide = AASequence::fromString(*sequence->second);
        }
        catch (const Exception::BaseException&)
        {
          require(false, "Could not parse library sequence: " + *sequence->second);
        }
        std::array<std::uint32_t, 3> range{};
        // Singly charged b_k and y_k for k = 1..n-1, kept inside the scan range used for the design study.
        for (int series = 0; series < 2; ++series)
        {
          range[series] = static_cast<std::uint32_t>(byion_mz_.size());
          for (Size k = 1; k < peptide.size(); ++k)
          {
            const double mz = series == 0 ? peptide.getPrefix(k).getMonoWeight(Residue::BIon, 1)
                                          : peptide.getSuffix(k).getMonoWeight(Residue::YIon, 1);
            if (mz >= 150.06 && mz <= 2000.0) { byion_mz_.push_back(static_cast<float>(mz)); }
          }
        }
        range[2] = static_cast<std::uint32_t>(byion_mz_.size());
        found = ranges.emplace(id, range).first;
      }
      byion_range_[i] = found->second;
    }
  }

  void TraceCollector::countByIons_(const MSSpectrum& spectrum, std::size_t index, double dt)
  {
    TraceEvidence& evidence = evidence_[index];
    // Spectra arrive in RT order, so the nearest isolating spectrum to the apex replaces any earlier one.
    if (!(dt < evidence.byion_apex_dt)) { return; }
    const auto& range = byion_range_[index];
    std::uint16_t counts[2] = {0, 0};
    for (std::uint32_t t = range[0]; t < range[2]; ++t)
    {
      const double mz = byion_mz_[t];
      const double half = mz * byion_ppm_ * 1e-6;
      double best = std::numeric_limits<double>::infinity();
      double best_intensity = 0.0;
      for (auto point = spectrum.MZBegin(mz - half); point != spectrum.end() && point->getMZ() <= mz + half; ++point)
      {
        const double distance = std::abs(point->getMZ() - mz);
        if (distance < best)
        {
          best = distance;
          best_intensity = point->getIntensity();
        }
      }
      if (std::isfinite(best) && best_intensity > byion_intensity_min_) { ++counts[t >= range[1] ? 1 : 0]; }
    }
    evidence.byion_apex_dt = static_cast<float>(dt);
    evidence.byion_b = counts[0];
    evidence.byion_y = counts[1];
  }

  void TraceCollector::consumeMS1_(MSSpectrum& spectrum)
  {
    const double rt = spectrum.getRT();
    auto first = std::lower_bound(by_rt_.begin(), by_rt_.end(), rt - rt_half_width_,
                                  [this](std::size_t i, double value) { return samples_[i].window.apex_rt < value; });
    if (first == by_rt_.end() || samples_[*first].window.apex_rt > rt + rt_half_width_)
    {
      return;
    }
    if (!spectrum.isSorted())
    {
      spectrum.sortByPosition();
    }
    // MS1 points are summed over ion mobility; the precursor m/z window is narrow enough for that on DIA data.
    const auto sum_window = [&](double target, double& ppm_sum, double& weight)
    {
      const double half = target * ms1_ppm_ * 1e-6;
      double total = 0.0;
      for (auto point = spectrum.MZBegin(target - half); point != spectrum.end() && point->getMZ() <= target + half; ++point)
      {
        const double intensity = point->getIntensity();
        if (!std::isfinite(intensity) || intensity <= 0.0) continue;
        total += intensity;
        ppm_sum += intensity * (point->getMZ() - target) / target * 1e6;
        weight += intensity;
      }
      return total;
    };
    std::vector<std::size_t> ms1_hits;
    for (auto active = first; active != by_rt_.end() && samples_[*active].window.apex_rt <= rt + rt_half_width_; ++active) ms1_hits.push_back(*active);
    // Each sample writes only its own evidence, so the samples are independent (thread-count invariant).
#pragma omp parallel for schedule(dynamic, 256) num_threads(threads_) if (threads_ > 1 && ms1_hits.size() >= 1024)
    for (std::ptrdiff_t h = 0; h < static_cast<std::ptrdiff_t>(ms1_hits.size()); ++h)
    {
      const std::size_t index = ms1_hits[static_cast<std::size_t>(h)];
      const Sample& sample = samples_[index];
      TraceEvidence& evidence = evidence_[index];
      const std::size_t bin = rtBin_(rt, sample.window.apex_rt);
      evidence.ms1_mask |= std::uint32_t{1} << bin;
      const bool core = std::abs(rt - sample.window.apex_rt) <= core_half_width_;
      double ppm_sum = 0.0, weight = 0.0;
      const double mono = sum_window(sample.precursor_mz, ppm_sum, weight);
      evidence.ms1[bin] += static_cast<float>(mono);
      if (core)
      {
        evidence.ms1_ppm_sum += ppm_sum;
        evidence.ms1_weight += weight;
        // First isotope (M+1) at the precursor charge, summed over the core window.
        if (sample.precursor_charge > 0)
        {
          double unused_sum = 0.0, unused_weight = 0.0;
          evidence.ms1_isotope += static_cast<float>(sum_window(sample.precursor_mz + 1.0033548 / sample.precursor_charge, unused_sum, unused_weight));
        }
      }
    }
  }

  void TraceCollector::consumeSpectrum(MSSpectrum& spectrum)
  {
    if (on_close_ && std::isfinite(spectrum.getRT()) && (spectrum.getMSLevel() == 1 || spectrum.getMSLevel() == 2)) close_(spectrum.getRT());
    if (spectrum.getMSLevel() == 1 && !evidence_.empty() && !spectrum.empty() && std::isfinite(spectrum.getRT()))
    {
      consumeMS1_(spectrum);
      return;
    }
    if (spectrum.getMSLevel() != 2 || spectrum.empty() || spectrum.getPrecursors().empty() ||
        !std::isfinite(spectrum.getRT()))
    {
      return;
    }
    const double rt = spectrum.getRT();
    auto first = std::lower_bound(by_rt_.begin(), by_rt_.end(), rt - rt_half_width_,
                                  [this](std::size_t i, double value) { return samples_[i].window.apex_rt < value; });
    if (first == by_rt_.end() || samples_[*first].window.apex_rt > rt + rt_half_width_)
    {
      return;
    }
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
    const bool scalar_valid = spectrum.getDriftTimeUnit() == DriftTimeUnit::VSSC &&
                               std::isfinite(scalar_im) && scalar_im > 0.0;
    // Isolation windows of this spectrum (same tests as per sample before; a NaN bound never isolates).
    std::vector<std::pair<double, double>> isolation;
    double lowest = std::numeric_limits<double>::infinity(), highest = -lowest;
    for (const auto& precursor : spectrum.getPrecursors())
    {
      const double lower = precursor.getIsolationWindowLowerOffset(), upper = precursor.getIsolationWindowUpperOffset();
      if (!(lower >= 0.0 && upper >= 0.0 && lower + upper > 0.0)) continue;
      const double from = precursor.getMZ() - lower, to = precursor.getMZ() + upper;
      if (std::isnan(from) || std::isnan(to)) continue;
      isolation.emplace_back(from, to);
      lowest = std::min(lowest, from);
      highest = std::max(highest, to);
    }
    if (isolation.empty()) return;
    const double last_bucket = static_cast<double>(by_mz_rt_.size() - 1);
    const double bucket_low = std::clamp(std::floor((lowest - mz_bucket_origin_) / MZ_BUCKET_TH), 0.0, last_bucket + 1.0);
    const double bucket_high = std::clamp(std::floor((highest - mz_bucket_origin_) / MZ_BUCKET_TH), -1.0, last_bucket);
    if (bucket_high < bucket_low) return;
    const std::size_t bucket_first = static_cast<std::size_t>(bucket_low);
    const std::size_t bucket_last = static_cast<std::size_t>(bucket_high);
    // Each sample sits in one bucket, so it is visited at most once per spectrum; per-sample state is independent of
    // the visiting order, so the result equals the former scan over all samples in by_rt_ order. The isolated samples
    // are collected (and their tensors assigned) serially, then filled in parallel: a sample's work writes only its
    // own tensor and evidence, in the same order as before, so the result does not depend on the thread count.
    hits_.clear();
    for (std::size_t bucket = bucket_first; bucket <= bucket_last; ++bucket)
    {
      const auto& in_bucket = by_mz_rt_[bucket];
      auto from_rt = std::lower_bound(in_bucket.begin(), in_bucket.end(), rt - rt_half_width_,
                                      [this](std::uint32_t i, double value) { return samples_[i].window.apex_rt < value; });
      for (auto active = from_rt; active != in_bucket.end() && samples_[*active].window.apex_rt <= rt + rt_half_width_; ++active)
      {
        const double precursor_mz = samples_[*active].precursor_mz;
        for (const auto& [from, to] : isolation)
        {
          if (from <= precursor_mz && to >= precursor_mz)
          {
            hits_.push_back(*active);
            hit_tensor_.resize(hits_.size());
            hit_tensor_.back() = tensor_(*active);
            break;
          }
        }
      }
    }
    std::exception_ptr error = nullptr;
#pragma omp parallel for schedule(dynamic, 64) num_threads(threads_) if (threads_ > 1 && hits_.size() >= 256)
    for (std::ptrdiff_t h = 0; h < static_cast<std::ptrdiff_t>(hits_.size()); ++h)
    {
    try
    {
      const std::size_t index = hits_[static_cast<std::size_t>(h)];
      Sample& sample = samples_[index];
      const bool compatible_im = sample.precursor_im > 0.0 && (mobilities != nullptr || scalar_valid);
      float* const tensor = hit_tensor_[static_cast<std::size_t>(h)];
      // Mixed availability cannot be represented by one sample mask. Collapse
      // such a sample to an XIC instead of inventing an IM position for missing data.
      if (!compatible_im && !sample.saw_missing_im && im_bins_ == IM_BINS)
      {
        for (std::size_t offset = 0; offset < SAMPLE_SIZE; offset += IM_BINS)
        {
          const double sum = std::accumulate(tensor + offset, tensor + offset + IM_BINS, 0.0);
          require(sum <= std::numeric_limits<float>::max(), "Collapsed trace intensity overflow.");
          std::fill_n(tensor + offset, IM_BINS, 0.0f);
          tensor[offset + IM_BINS / 2] = static_cast<float>(sum);
        }
        sample.saw_missing_im = true;
        sample.has_im = false;
      }
      if (compatible_im && !sample.saw_missing_im)
      {
        sample.has_im = true;
      }
      const auto rt_bin = std::min(RT_BINS - 1, static_cast<std::size_t>(
        std::max(0.0, (rt - sample.window.apex_rt + rt_half_width_) / (2.0 * rt_half_width_)) * RT_BINS));
      TraceEvidence* evidence = evidence_.empty() ? nullptr : &evidence_[index];
      const bool core = evidence != nullptr && std::abs(rt - sample.window.apex_rt) <= core_half_width_;
      if (evidence != nullptr)
      {
        evidence->ms2_mask |= std::uint32_t{1} << rt_bin;
        if (core && !byion_range_.empty()) { countByIons_(spectrum, index, std::abs(rt - sample.window.apex_rt)); }
      }
      for (std::size_t fragment = 0; fragment < FRAGMENTS; ++fragment)
      {
        const double mz = sample.fragment_mz[fragment];
        if (mz == 0.0)
        {
          continue;
        }
        for (auto point = spectrum.MZBegin(mz - mz_tolerance_);
             point != spectrum.end() && point->getMZ() <= mz + mz_tolerance_; ++point)
        {
          const double intensity = point->getIntensity();
          if (!std::isfinite(intensity) || intensity <= 0.0)
          {
            continue;
          }
          std::size_t im_bin = IM_BINS / 2;
          if (compatible_im)
          {
            const double im = mobilities != nullptr ? (*mobilities)[point - spectrum.begin()] : scalar_im;
            if (!std::isfinite(im) || im <= 0.0 || std::abs(im - sample.precursor_im) > im_half_width_)
            {
              continue;
            }
            if (!sample.saw_missing_im)
            {
              im_bin = std::min(IM_BINS - 1, static_cast<std::size_t>(
                std::max(0.0, (im - sample.precursor_im + im_half_width_) / (2.0 * im_half_width_)) * IM_BINS));
            }
          }
          float& value = im_bins_ == IM_BINS ? tensor[(fragment * RT_BINS + rt_bin) * IM_BINS + im_bin]
                                             : tensor[fragment * RT_BINS + rt_bin];
          require(intensity <= std::numeric_limits<float>::max() - static_cast<double>(value), "Trace intensity overflow.");
          value += static_cast<float>(intensity);
          if (core)
          {
            evidence->ms2_ppm_sum += intensity * (point->getMZ() - mz) / mz * 1e6;
            evidence->ms2_weight += intensity;
          }
        }
        if (core && !pg_evidence_.empty())
        {
          // Sum of the points within pg_ppm_ of target (and the tensor's mobility window); optionally their ppm sum.
          const auto sum_at = [&](double target, double* ppm_sum)
          {
            const double tol = target * pg_ppm_ * 1e-6;
            double total = 0.0;
            for (auto point = spectrum.MZBegin(target - tol); point != spectrum.end() && point->getMZ() <= target + tol; ++point)
            {
              const double intensity = point->getIntensity();
              if (!std::isfinite(intensity) || intensity <= 0.0) continue;
              if (compatible_im)
              {
                const double im = mobilities != nullptr ? (*mobilities)[point - spectrum.begin()] : scalar_im;
                if (!std::isfinite(im) || im <= 0.0 || std::abs(im - sample.precursor_im) > im_half_width_) continue;
              }
              total += intensity;
              if (ppm_sum != nullptr) *ppm_sum += intensity * (point->getMZ() - target) / target * 1e6;
            }
            return total;
          };
          auto& pg = pg_evidence_[index];
          double ppm_sum = 0.0;
          pg.mono[fragment] += static_cast<float>(sum_at(mz, &ppm_sum));
          pg.ppm_sum[fragment] += static_cast<float>(ppm_sum);
          const double z = std::abs(sample.fragment_charge[fragment] != 0 ? static_cast<double>(sample.fragment_charge[fragment]) : 1.0);
          for (std::size_t iso = 0; iso < PeakGroupEvidence::ISOTOPES; ++iso)
          {
            pg.isotope[fragment][iso] += static_cast<float>(sum_at(mz + (iso + 1) * 1.0033548 / z, nullptr));
            pg.before[fragment][iso] += static_cast<float>(sum_at(mz - 1.0033548 / static_cast<double>(iso + 1), nullptr));
          }
        }
      }
    }
    catch (...)
    {
#pragma omp critical (OpenDIA_TraceCollector_error)
      {
        if (!error) error = std::current_exception();
      }
    }
    }
    if (error) std::rethrow_exception(error);
  }

  float* TraceCollector::tensor_(std::size_t index)
  {
    if (!on_close_) return intensities_.data() + index * stride_;
    require(!closed_[index], "TraceCollector: spectrum out of RT order reaches a closed window.");
    std::uint32_t& slot = slot_of_[index];
    if (slot == NO_SLOT)
    {
      require(held_ < max_held_, "Trace windows held at one time exceed the raw tensor byte budget.");
      if (free_slots_.empty())
      {
        const std::size_t base = chunks_.size() * SLOT_CHUNK;
        chunks_.emplace_back(new float[SLOT_CHUNK * stride_]());
        for (std::size_t s = SLOT_CHUNK; s-- > 0;) free_slots_.push_back(static_cast<std::uint32_t>(base + s));
      }
      slot = free_slots_.back();
      free_slots_.pop_back();
      peak_held_ = std::max(peak_held_, ++held_);
    }
    return chunks_[slot / SLOT_CHUNK].get() + (slot % SLOT_CHUNK) * stride_;
  }

  const float* TraceCollector::storedTensor_(std::size_t index) const
  {
    if (!on_close_) return intensities_.data() + index * stride_;
    const std::uint32_t slot = slot_of_[index];
    return slot == NO_SLOT ? nullptr : chunks_[slot / SLOT_CHUNK].get() + (slot % SLOT_CHUNK) * stride_;
  }

  void TraceCollector::close_(double rt)
  {
    // A window receives spectra with RT in [apex - half, apex + half]; once rt passes apex + half it is complete.
    while (close_cursor_ < by_rt_.size() && samples_[by_rt_[close_cursor_]].window.apex_rt < rt - rt_half_width_)
    {
      const auto index = static_cast<std::uint32_t>(by_rt_[close_cursor_++]);
      closed_[index] = 1;
      pending_.push_back(index);
    }
    if (pending_.size() >= CLOSE_BATCH) flush_();
  }

  void TraceCollector::flush_()
  {
    if (pending_.empty()) return;
    on_close_(pending_, *this);
    for (const std::uint32_t index : pending_)
    {
      std::uint32_t& slot = slot_of_[index];
      if (slot == NO_SLOT) continue;
      std::fill_n(chunks_[slot / SLOT_CHUNK].get() + (slot % SLOT_CHUNK) * stride_, stride_, 0.0f);
      free_slots_.push_back(slot);
      slot = NO_SLOT;
      --held_;
    }
    pending_.clear();
  }

  void TraceCollector::finish()
  {
    if (!on_close_) return;
    close_(std::numeric_limits<double>::infinity());
    flush_();
    std::vector<std::unique_ptr<float[]>>().swap(chunks_);
    std::vector<std::uint32_t>().swap(free_slots_);
  }

  const float* TraceCollector::windowTensor(std::size_t i, std::array<float, SAMPLE_SIZE>& buffer) const
  {
    require(i < samples_.size(), "Window index out of range.");
    const float* stored = storedTensor_(i);
    if (stored == nullptr)
    {
      buffer.fill(0.0f);
      return buffer.data();
    }
    if (im_bins_ == IM_BINS) { return stored; }
    buffer.fill(0.0f);
    for (std::size_t k = 0; k < FRAGMENTS * RT_BINS; ++k) { buffer[k * IM_BINS + IM_BINS / 2] = stored[k]; }
    return buffer.data();
  }

  void TraceCollector::consumeChromatogram(MSChromatogram&) {}
  void TraceCollector::setExpectedSize(std::size_t, std::size_t) {}
  void TraceCollector::setExperimentalSettings(const ExperimentalSettings&) {}
}
