// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace OpenMS::OpenDIA
{
  /**
    @brief Shared-fragment explaining away between co-eluting candidate windows (the interference graph, one pass).

    A window is one candidate peak of one precursor. Per window the resolver needs the peak range [rt_start, rt_end],
    the apex, the central area of each fragment and the summed fragment profile on a fixed RT grid of BINS bins over
    apex +- half width (the grid of the 12-fragment features). Per precursor it needs its isolation window and up to
    FRAGMENTS fragment m/z with library intensities.

    Two windows a and b of different precursors are linked when their precursors share an isolation window, their peak
    ranges overlap, at least one fragment of a lies within the fragment tolerance (ppm of a's fragment m/z) of one
    of b, and (with im_link_width > 0 and both 1/K0 known) their 1/K0 differ by at most im_link_width. For each
    (a, claimant precursor) only the claimant window with the nearest apex is used. The claimant b
    explains part of a's signal on each shared fragment f (matched to b's fragment g):

      scale_b = median over b's fragments not shared with a, with library intensity and central area > 0, of area / lib
      frac_ab = b's summed profile over a's central bins (shifted by the apex difference in bins, rounded half to even)
                / b's summed profile over its own central bins
      expect  = scale_b * lib_b[g] * frac_ab

    No unshared fragment with signal means b explains nothing. Expectations of all claimants are summed per fragment of a
    and capped at a's observed area. Everything is identical for targets and decoys and uses no labels.

    Features per window, in featureNames() order:
      i_explained     explained share of the window's central area (0 without signal)
      i_libcos_resid  cosine of sqrt(residual area) and sqrt(library intensity) over the fragments
      i_dlibcos       i_libcos_resid minus the same cosine on the observed areas
      i_nclaim        log(1 + claimant precursors with a nonzero contribution)

    Reproduces gap-20261114/scripts/features_I.py (PLAN_GAP arm I, validated on the Astral entrapment runs).
  */
  class InterferenceResolver
  {
  public:
    static constexpr std::size_t FRAGMENTS = 12;
    static constexpr std::size_t BINS = 32;
    /// Central bins [CORE_BEGIN, CORE_END) of the grid: apex +- 4.7 s on the +-15 s grid of the 12-fragment features.
    static constexpr std::size_t CORE_BEGIN = 12;
    static constexpr std::size_t CORE_END = 20;
    static constexpr std::size_t FEATURES = 4;

    struct Precursor
    {
      std::int32_t isolation = -1;                 ///< Isolation window index; -1 never links.
      std::array<double, FRAGMENTS> fragment_mz{};  ///< 0 marks an absent fragment.
      std::array<double, FRAGMENTS> library{};      ///< Library intensity relative to the precursor's most intense fragment; 0 for absent fragments.
      double im = 0.0;                              ///< Calibrated 1/K0; 0 = unknown (links on m/z only).
    };

    struct Window
    {
      std::uint32_t precursor = 0;                  ///< Index into the precursor array.
      double rt_start = 0.0;
      double rt_end = 0.0;
      double rt_apex = 0.0;
      std::array<double, FRAGMENTS> core_area{};    ///< Central area per fragment slot.
      std::array<double, BINS> profile{};           ///< Summed fragment intensity per grid bin.
    };

    struct Statistics
    {
      std::size_t pairs_tested = 0;          ///< Overlapping window pairs of different precursors in one isolation window.
      std::size_t links = 0;                 ///< Directed links with at least one shared fragment.
      std::size_t used_links = 0;            ///< After keeping the nearest claimant window per precursor.
      std::size_t windows_with_claimant = 0;
    };

    /**
      @brief Compute the features of every window.
      @param[in] precursors Precursors referenced by the windows.
      @param[in] windows Candidate windows.
      @param[in] fragment_tolerance_ppm Fragment match half width in ppm (> 0).
      @param[in] rt_bin_seconds Width of one grid bin in seconds (> 0).
      @param[in] threads Worker threads (>= 1); the result does not depend on it.
      @param[out] stats Optional counts.
      @param[in] im_link_width If > 0, two precursors that both have a 1/K0 link only if their 1/K0 differ by at most this
                 (their mobility extraction windows overlap); 0 links on m/z alone.
      @return windows.size() x FEATURES, row-major.
      @throws std::invalid_argument On a non-positive tolerance, bin width or thread count, or a bad precursor index.
    */
    static std::vector<float> compute(const std::vector<Precursor>& precursors, const std::vector<Window>& windows,
                                      double fragment_tolerance_ppm, double rt_bin_seconds, int threads = 1, Statistics* stats = nullptr,
                                      double im_link_width = 0.0);

    /// Feature column names, in output order.
    static const std::array<const char*, FEATURES>& featureNames();
  };
}
