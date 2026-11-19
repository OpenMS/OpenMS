// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenDIAInterference.h>

#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

namespace OpenMS::OpenDIA
{
  /**
    @brief Up to 12 library fragments of each candidate window on a 32-bin RT grid (apex +- half width), and the
    label-free evidence computed from it.

    A grid holds, per fragment slot f and bin b, the summed intensity of every MS2 point within the extraction ppm of the
    fragment's calibrated m/z, over the spectra that isolate the precursor. The central bins [12, 20) are apex +- 4.7 s on
    the +-15 s grid. Slots follow the scoring library order; `library` is the intensity relative to the precursor's most
    intense fragment (0 = no fragment) and `search` marks fragments that are also in the search library.
  */
  struct FragmentGrid
  {
    static constexpr std::size_t FRAGMENTS = InterferenceResolver::FRAGMENTS;
    static constexpr std::size_t BINS = InterferenceResolver::BINS;
    static constexpr std::size_t CORE_BEGIN = InterferenceResolver::CORE_BEGIN;
    static constexpr std::size_t CORE_END = InterferenceResolver::CORE_END;
    static constexpr std::size_t EVIDENCE = 6;
    static constexpr std::size_t DEFICIT = 3;
    using Values = std::array<double, FRAGMENTS>;
    using Flags = std::array<bool, FRAGMENTS>;

    /// Column names of evidence() and deficit(), in output order.
    static const std::array<const char*, EVIDENCE>& evidenceNames();
    static const std::array<const char*, DEFICIT>& deficitNames();

    /**
      @brief Six fragment-evidence features of one grid ([FRAGMENTS][BINS], row-major):

      f12_support          fragments with central signal whose profile correlates > 0.5 with the sum of the others
      f12_coelution        mean of those correlations over the library fragments
      f12_libcos           cosine of sqrt(central area) and sqrt(library intensity)
      f12_extra_support    f12_support over fragments not in the search library
      f12_extra_coelution  mean correlation over those fragments
      f12_extra_share      log10((their central area + 1) / (all central area + 1))
    */
    static std::array<float, EVIDENCE> evidence(const double* grid, const Values& library, const Flags& search);

    /**
      @brief Fragment-deficit features d_max, d6_max, d_lib of one window.

      With scale the median of area_f / library_f over library fragments with central signal, the deficit of fragment f
      is log10(scale library_f + noise) - log10(area_f + noise): how far it falls below the abundance its siblings imply.
      d_max is the largest deficit over the library fragments, d6_max over the search-library fragments, and d_lib the
      library-intensity share of fragments with a deficit above 1.
      @param[in] noise The run's noise floor (noiseLevel()).
    */
    static std::array<float, DEFICIT> deficit(const Values& core_area, const Values& library, const Flags& search, double noise);

    /// Median of the nonzero central areas of library fragments over all windows (1 if there are none).
    static double noiseLevel(const std::vector<InterferenceResolver::Window>& windows, const std::vector<InterferenceResolver::Precursor>& precursors);
  };

  /**
    @brief Stream MS2 spectra once and fill the fragment grid of every candidate window.

    For every window and every MS2 spectrum isolating the window's precursor with RT inside apex +- rt_half_width
    (open interval), the intensity of every point within extraction_ppm of each calibrated fragment m/z goes to bin
    floor((rt - apex + half) / (2 half) x BINS). A scan within 1 ms of the apex goes to the upper central bin. Ion
    mobility is ignored (points are pooled). Fragment sums are computed once per (precursor, spectrum).

    Grids are held only while a window is open. When the spectrum RT passes apex + half the window is closed: its central
    areas and summed profile are stored in the Window (the interference graph inputs) and its evidence features are
    computed. Spectra must arrive in RT order; an earlier spectrum after a later one throws.

    Each precursor is assigned the distinct isolation window with the largest lower edge at or below its m/z, if its m/z
    lies below that window's upper edge.
  */
  class FragmentGridCollector final : public Interfaces::IMSDataConsumer
  {
  public:
    struct Precursor
    {
      double mz = 0.0;
      std::array<double, FragmentGrid::FRAGMENTS> extraction_mz{}; ///< Calibrated fragment m/z per slot (0 = absent).
      FragmentGrid::Values library{};
      FragmentGrid::Flags search{};
      double im = 0.0; ///< Calibrated 1/K0; 0 = no mobility filter for this precursor.
    };

    /**
      @param[in] precursors Indexed by InterferenceResolver::Window::precursor.
      @param[in,out] windows rt_apex set by the caller; core_area and profile are filled.
      @param[in] rt_half_width Half width of the grid in seconds (> 0).
      @param[in] extraction_ppm Half width of each fragment extraction in ppm (> 0).
    */
    FragmentGridCollector(const std::vector<Precursor>& precursors, std::vector<InterferenceResolver::Window>& windows, double rt_half_width,
                          double extraction_ppm);

    /**
      @brief Restrict extraction to points within @p half_width (1/K0) of each precursor's calibrated 1/K0.

      Only spectra with a VSSC mobility array (or a scalar VSSC drift time), and only precursors with im > 0, are
      filtered; everything else keeps pooled extraction. Must be called before the first spectrum.
      @param[in] half_width Half width in 1/K0 (> 0).
      @param[in] array_as_vssc Treat a mobility array whose unit is not CCS as 1/K0.
    */
    void enableIonMobility(double half_width, bool array_as_vssc = false);

    void consumeSpectrum(MSSpectrum& spectrum) override;
    void consumeChromatogram(MSChromatogram&) override {}
    void setExpectedSize(std::size_t, std::size_t) override {}
    void setExperimentalSettings(const ExperimentalSettings&) override {}

    /// Close every window still open; call once after the last spectrum.
    void finish();

    /// Evidence features per window (FragmentGrid::evidence()), valid after finish().
    const std::vector<std::array<float, FragmentGrid::EVIDENCE>>& evidence() const { return evidence_; }
    /// Isolation window index per precursor (-1 none).
    std::vector<std::int32_t> isolationIndices() const;
    std::size_t isolationWindowCount() const { return lanes_.size(); }
    std::size_t spectraUsed() const { return spectra_used_; }
    /// Largest number of grids held at once.
    std::size_t peakOpenGrids() const { return peak_open_; }
    /// Worker threads for the per-spectrum fragment sums (results do not depend on it). Default 1.
    void setThreads(int threads) { threads_ = threads < 1 ? 1 : threads; }

  private:
    struct Lane
    {
      std::vector<std::uint32_t> precursors; ///< Precursors with m/z in the isolation window, by first window apex - half.
      std::size_t cursor = 0;
      std::vector<std::uint32_t> active;
    };

    Lane& lane_(double lo, double hi);
    double* grid_(std::uint32_t window);
    void close_(std::uint32_t window);
    void closeUntil_(double rt);

    const std::vector<Precursor>& precursors_;
    std::vector<InterferenceResolver::Window>& windows_;
    std::vector<std::uint32_t> window_offsets_, window_index_;
    std::vector<double> rt_low_, rt_high_;
    std::vector<std::uint32_t> by_apex_;  ///< Windows by apex (close order).
    std::size_t close_cursor_ = 0;
    std::vector<std::int32_t> slot_;      ///< Grid slot per window (-1 none).
    std::vector<double> pool_;
    std::vector<std::int32_t> free_;
    std::vector<std::array<float, FragmentGrid::EVIDENCE>> evidence_;
    std::map<std::pair<double, double>, Lane> lanes_;
    double half_, ppm_;
    double im_half_ = 0.0;
    bool im_array_as_vssc_ = false;
    double last_rt_ = -1e300;
    std::size_t spectra_used_ = 0, open_ = 0, peak_open_ = 0;
    bool finished_ = false;
    int threads_ = 1;
    std::vector<std::uint32_t> reached_;                         ///< Precursors reached by the current spectrum.
    std::vector<std::array<double, FragmentGrid::FRAGMENTS>> sums_; ///< Their fragment sums.
    std::vector<std::uint8_t> any_;
  };
}
