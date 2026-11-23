// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace OpenMS::OpenDIA
{
  /**
    @brief Piecewise-linear library-RT to run-RT map from signal anchors.

    Knots are the median (library RT, run RT) of 20 equal-count bins over the
    anchors sorted by library RT. Between knots the map is piecewise linear.
    Outside the knot range it extrapolates linearly, on the line through the end
    knot and the knot kEndSpan positions inward (slope clamped at zero, so the map
    stays nondecreasing). Held-out end anchors on the Astral entrapment runs: the
    median absolute residual is 17-22 s with linear ends, against 67-69 s with
    constant ends.
  */
  struct RTMap
  {
    /// Knot offset from each end that sets the extrapolation slope.
    static constexpr std::size_t kEndSpan = 3;
    std::vector<double> library_rt; ///< Knot library RT coordinates, nondecreasing.
    std::vector<double> run_rt; ///< Knot run RT in seconds.
    double sigma = 0.0; ///< 1.4826 x median absolute residual of the final fit, in seconds.
    std::size_t anchors = 0; ///< Anchors kept after outlier rejection.
    /// Optional extrapolation slopes (s per library RT unit) below the first and above the last knot. NaN uses the
    /// kEndSpan rule.
    double low_slope = std::numeric_limits<double>::quiet_NaN();
    double high_slope = std::numeric_limits<double>::quiet_NaN();

    /// Predicted run RT in seconds for a library RT coordinate.
    double operator()(double library_rt) const;
  };

  /**
    @brief Fit an RTMap with up to three outlier-rejection refits.

    Each refit drops anchors whose absolute residual exceeds max(3 sigma, 3 s).

    @param[in] anchors (library RT, run RT in seconds) pairs; at least 20 finite pairs.
    @param[in] bins Equal-count bins (knots); default 20.
    @return The fitted map.
    @throws std::invalid_argument If there are too few or non-finite anchors.
  */
  RTMap fitRTMap(std::vector<std::pair<double, double>> anchors, std::size_t bins = 20);

  /**
    @brief Fit a nondecreasing RTMap without a bin count: L1 isotonic regression of run RT on library RT.

    Each block of the isotonic fit becomes one knot at (median library RT, block median run RT), so the knot
    density follows the anchors and isolated anchors at the gradient ends keep their own knots instead of being
    averaged into a wide end bin. Beyond the knots the map extrapolates with the Theil-Sen slope of the
    round(sqrt(n)) anchors at that end (at least 3; clamped at zero). sigma is 1.4826 x the median absolute
    10-fold cross-validated residual (interleaved folds over the sorted anchors), the prediction error of a
    precursor not used in the fit. Up to three refits drop anchors whose cross-validated residual exceeds
    max(3 sigma, 3 s).

    @param[in] anchors (library RT, run RT in seconds) pairs; at least 20 finite pairs.
    @return The fitted map.
    @throws std::invalid_argument If there are too few or non-finite anchors, or rejection leaves fewer than 20.
  */
  RTMap fitMonotoneRTMap(std::vector<std::pair<double, double>> anchors);
}
