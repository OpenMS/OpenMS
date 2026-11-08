// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIAPeakMaps.h"

#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/METADATA/Precursor.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace OpenMS::OpenDIA
{
  namespace
  {
    constexpr double C13_C12_MASS_DIFFERENCE = 1.0033548;

    /// Per-peak or scalar 1/K0 of a spectrum; the same rules as StreamSearch.
    struct SpectrumMobility
    {
      const DataArrays::FloatDataArray* per_peak = nullptr;
      double scalar = -1.0;
    };

    SpectrumMobility spectrumMobility(const MSSpectrum& spectrum, bool array_as_vssc)
    {
      SpectrumMobility mobility;
      if (spectrum.containsIMData())
      {
        const auto [im_index, im_unit] = spectrum.getIMData();
        const auto& arrays = spectrum.getFloatDataArrays();
        if (im_index < arrays.size() && arrays[im_index].size() == spectrum.size() &&
            (im_unit == DriftTimeUnit::VSSC || (array_as_vssc && im_unit != DriftTimeUnit::CCS)))
        {
          mobility.per_peak = &arrays[im_index];
        }
      }
      if (mobility.per_peak == nullptr && spectrum.getDriftTime() != IMTypes::DRIFTTIME_NOT_SET &&
          spectrum.getDriftTimeUnit() == DriftTimeUnit::VSSC && std::isfinite(spectrum.getDriftTime()) &&
          spectrum.getDriftTime() > 0.0)
      {
        mobility.scalar = spectrum.getDriftTime();
      }
      return mobility;
    }

    bool isolationContains(const MSSpectrum& spectrum, double mz)
    {
      for (const Precursor& precursor : spectrum.getPrecursors())
      {
        const double center = precursor.getMZ();
        const double lower = precursor.getIsolationWindowLowerOffset();
        const double upper = precursor.getIsolationWindowUpperOffset();
        if (std::isfinite(center) && std::isfinite(lower) && std::isfinite(upper) && lower >= 0.0 && upper >= 0.0 &&
            lower + upper > 0.0 && mz >= center - lower && mz <= center + upper)
        {
          return true;
        }
      }
      return false;
    }
  }

  PeakMapCollector::PeakMapCollector(std::vector<PeakMapTarget> targets, PeakMapOptions options, Sink sink) :
    targets_(std::move(targets)), options_(options), sink_(std::move(sink))
  {
    if (!sink_ || !std::isfinite(options_.mz_tolerance_ppm) || !(options_.mz_tolerance_ppm > 0.0) ||
        !std::isfinite(options_.rt_half_window_s) || options_.rt_half_window_s < 0.0 ||
        !std::isfinite(options_.im_half_window) || options_.im_half_window < 0.0)
    {
      throw std::invalid_argument("Peak-map export needs a sink, a positive m/z tolerance and non-negative RT and IM windows.");
    }
    states_.resize(targets_.size());
    for (std::size_t i = 0; i < targets_.size(); ++i)
    {
      const PeakMapTarget& target = targets_[i];
      State& state = states_[i];
      const bool windowed = target.has_rt && std::isfinite(target.center_rt) && options_.rt_half_window_s > 0.0;
      if (windowed)
      {
        state.rt_start = target.center_rt - options_.rt_half_window_s;
        state.rt_end = target.center_rt + options_.rt_half_window_s;
      }
      const double target_im = target.has_im ? target.precursor_im : -1.0;
      // MS1: monoisotope plus isotopes. Without a charge only the monoisotope is defined.
      const std::size_t isotopes = target.charge > 0 ? options_.ms1_isotopes : 0;
      for (std::size_t k = 0; k <= isotopes; ++k)
      {
        PeakMapExtractor::ExtractedPeakMap map;
        map.native_id = target.precursor_id + "_Precursor_i" + std::to_string(k);
        map.target_mz = target.precursor_mz + (target.charge > 0 ? k * C13_C12_MASS_DIFFERENCE / target.charge : 0.0);
        map.target_rt = target.has_rt ? target.center_rt : std::numeric_limits<double>::quiet_NaN();
        map.target_ion_mobility = target_im;
        map.rt_start = windowed ? state.rt_start : 0.0;
        map.rt_end = windowed ? state.rt_end : -1.0;
        state.ms1.push_back(std::move(map));
      }
      for (const auto& [native_id, fragment_mz] : target.fragments)
      {
        PeakMapExtractor::ExtractedPeakMap map;
        map.native_id = native_id;
        map.target_mz = fragment_mz;
        map.target_rt = target.has_rt ? target.center_rt : std::numeric_limits<double>::quiet_NaN();
        map.target_ion_mobility = target_im;
        map.rt_start = windowed ? state.rt_start : 0.0;
        map.rt_end = windowed ? state.rt_end : -1.0;
        state.ms2.push_back(std::move(map));
      }
    }
    open_order_.resize(targets_.size());
    std::iota(open_order_.begin(), open_order_.end(), std::size_t{0});
    std::stable_sort(open_order_.begin(), open_order_.end(),
                     [this](std::size_t lhs, std::size_t rhs) { return states_[lhs].rt_start < states_[rhs].rt_start; });
  }

  void PeakMapCollector::consumeSpectrum(MSSpectrum& spectrum)
  {
    if (finished_ || spectrum.empty())
    {
      return;
    }
    const int level = static_cast<int>(spectrum.getMSLevel());
    if (level != 1 && (level != 2 || spectrum.getPrecursors().empty()))
    {
      return;
    }
    const double rt = spectrum.getRT();
    if (!std::isfinite(rt))
    {
      return;
    }
    if (rt < max_rt_)
    {
      ++out_of_order_;
    }
    else
    {
      max_rt_ = rt;
      // Open targets whose window has started; close targets the run has passed.
      while (next_open_ < open_order_.size() && states_[open_order_[next_open_]].rt_start <= rt)
      {
        active_.push_back(open_order_[next_open_++]);
      }
      std::size_t kept = 0;
      for (std::size_t target : active_)
      {
        if (states_[target].rt_end < rt)
        {
          emit_(target);
        }
        else
        {
          active_[kept++] = target;
        }
      }
      active_.resize(kept);
    }
    extract_(spectrum, rt);
  }

  void PeakMapCollector::extract_(const MSSpectrum& spectrum, double rt)
  {
    const bool ms1 = spectrum.getMSLevel() == 1;
    const SpectrumMobility mobility = spectrumMobility(spectrum, options_.ion_mobility_array_as_vssc);
    const bool sorted = std::is_sorted(spectrum.begin(), spectrum.end(),
                                       [](const Peak1D& lhs, const Peak1D& rhs) { return lhs.getMZ() < rhs.getMZ(); });
    const auto collect = [&](PeakMapExtractor::ExtractedPeakMap& map, const PeakMapTarget& target)
    {
      const double half = map.target_mz * options_.mz_tolerance_ppm * 1e-6;
      const double low = map.target_mz - half;
      const double high = map.target_mz + half;
      std::size_t begin = 0;
      if (sorted)
      {
        begin = static_cast<std::size_t>(std::lower_bound(spectrum.begin(), spectrum.end(), low,
                                                          [](const Peak1D& peak, double mz) { return peak.getMZ() < mz; }) -
                                         spectrum.begin());
      }
      for (std::size_t i = begin; i < spectrum.size(); ++i)
      {
        const double mz = spectrum[i].getMZ();
        if (mz > high)
        {
          if (sorted) break;
          continue;
        }
        if (mz < low)
        {
          continue;
        }
        double im = mobility.per_peak != nullptr ? static_cast<double>((*mobility.per_peak)[i]) : mobility.scalar;
        if (!std::isfinite(im) || im <= 0.0)
        {
          im = -1.0;
        }
        if (options_.im_half_window > 0.0 && target.has_im && im > 0.0 &&
            std::abs(im - target.precursor_im) > options_.im_half_window)
        {
          continue;
        }
        map.mz.push_back(mz);
        map.rt.push_back(rt);
        map.ion_mobility.push_back(im);
        map.intensity.push_back(spectrum[i].getIntensity());
      }
    };
    for (std::size_t target_index : active_)
    {
      State& state = states_[target_index];
      if (rt < state.rt_start || rt > state.rt_end || state.emitted)
      {
        continue;
      }
      const PeakMapTarget& target = targets_[target_index];
      if (ms1)
      {
        for (auto& map : state.ms1) collect(map, target);
      }
      else if (isolationContains(spectrum, target.precursor_mz))
      {
        for (auto& map : state.ms2) collect(map, target);
      }
    }
  }

  void PeakMapCollector::emit_(std::size_t target)
  {
    State& state = states_[target];
    if (state.emitted)
    {
      return;
    }
    state.emitted = true;
    for (const auto& map : state.ms1)
    {
      sink_(map, 1);
      points_written_ += map.mz.size();
      ++peak_maps_written_;
    }
    for (const auto& map : state.ms2)
    {
      sink_(map, 2);
      points_written_ += map.mz.size();
      ++peak_maps_written_;
    }
    // Release the points; the emitted flag keeps late spectra out.
    state.ms1 = {};
    state.ms2 = {};
  }

  void PeakMapCollector::finish()
  {
    if (finished_)
    {
      return;
    }
    finished_ = true;
    // Active targets first, then targets whose window never started, so every target is written once.
    for (std::size_t target : active_) emit_(target);
    active_.clear();
    for (; next_open_ < open_order_.size(); ++next_open_) emit_(open_order_[next_open_]);
  }

  void PeakMapCollector::consumeChromatogram(MSChromatogram&) {}
  void PeakMapCollector::setExpectedSize(std::size_t, std::size_t) {}
  void PeakMapCollector::setExperimentalSettings(const ExperimentalSettings&) {}
  std::size_t PeakMapCollector::outOfOrderSpectra() const { return out_of_order_; }
  std::size_t PeakMapCollector::pointsWritten() const { return points_written_; }
  std::size_t PeakMapCollector::peakMapsWritten() const { return peak_maps_written_; }

  TeeConsumer::TeeConsumer(Interfaces::IMSDataConsumer& first, Interfaces::IMSDataConsumer& second) :
    first_(first), second_(second)
  {
  }

  void TeeConsumer::consumeSpectrum(MSSpectrum& spectrum)
  {
    first_.consumeSpectrum(spectrum);
    second_.consumeSpectrum(spectrum);
  }

  void TeeConsumer::consumeChromatogram(MSChromatogram& chromatogram)
  {
    first_.consumeChromatogram(chromatogram);
    second_.consumeChromatogram(chromatogram);
  }

  void TeeConsumer::setExpectedSize(std::size_t spectra, std::size_t chromatograms)
  {
    first_.setExpectedSize(spectra, chromatograms);
    second_.setExpectedSize(spectra, chromatograms);
  }

  void TeeConsumer::setExperimentalSettings(const ExperimentalSettings& settings)
  {
    first_.setExperimentalSettings(settings);
    second_.setExperimentalSettings(settings);
  }
}
