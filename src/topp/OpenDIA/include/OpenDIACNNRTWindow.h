// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenDIACNNQValues.h>
#include <OpenDIARTMap.h>

#include <cstddef>
#include <utility>
#include <vector>

namespace OpenMS::OpenDIA
{
  /** @brief Per-fold summary of the out-of-fold RT window. */
  struct RTWindowFold
  {
    int fold = 0;
    RTMap map; ///< Fitted only on anchors from the other folds.
    double half_window = 0.0; ///< max(4 sigma, 30 s), in seconds.
    std::size_t precursors = 0; ///< Precursors in this fold.
    std::size_t kept_precursors = 0; ///< Precursors with at least one window inside the RT window.
  };

  /**
    @brief Experimental signal-derived RT window for crossfit scoring.

    Anchors are the @p anchors_per_fold target precursors with the highest best-window
    score in each fold (ties: larger precursor ID), at their best window's apex. No
    teacher identifications are used. The map for fold k is fitted only on anchors
    from the other folds. A window is kept iff |apex - map(library RT)| <= max(4 sigma, 30 s).

    @param[in] inputs Out-of-fold scored windows (pair keys are ignored).
    @param[in] apex_rt Run RT in seconds of each input window's apex.
    @param[in] library_rt Library RT coordinate of each input window's precursor.
    @param[out] report One summary per represented fold, in fold order.
    @param[in] anchors_per_fold Anchor targets taken from each fold.
    @return Mask over @p inputs; true keeps the window.
    @throws std::invalid_argument On size mismatch, non-finite coordinates, fewer than two folds or too few anchors.
  */
  std::vector<bool> crossfitRTWindowMask(const std::vector<CNNQValueInput>& inputs, const std::vector<double>& apex_rt,
                                         const std::vector<double>& library_rt, std::vector<RTWindowFold>& report,
                                         std::size_t anchors_per_fold = 500);
}
