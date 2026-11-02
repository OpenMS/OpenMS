// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <algorithm>
#include <array>
#include <limits>
#include <string_view>
#include <vector>

namespace OpenMS
{
  /// Helpers of NuXL:rescoring: the delta score, the training hits of the self-trained ion likelihood, the native score
  struct NuXLRescoring
  {
    /**
      @brief Best score of the hits of a spectrum and the best score of a different (unmodified) peptide sequence.

      Fed every hit scored for the spectrum in the search, in any order. bestOtherThan() then gives the exact best
      score of a sequence other than the given one, independent of the order of the hits.
    */
    struct ScoreGap
    {
      float best = -std::numeric_limits<float>::infinity();
      float second = -std::numeric_limits<float>::infinity(); ///< best score of a sequence other than best_sequence
      std::string_view best_sequence;

      void add(float score, std::string_view sequence)
      {
        if (sequence == best_sequence) { best = std::max(best, score); return; }
        if (score > best)
        {
          second = best; // the old best has another sequence and is at least as good as every other hit seen so far
          best = score;
          best_sequence = sequence;
          return;
        }
        second = std::max(second, score);
      }

      /// best score of a sequence other than @p sequence (-inf if there is none)
      float bestOtherThan(std::string_view sequence) const { return sequence == best_sequence ? second : best; }
    };

    /// FDR threshold of the training hits of the self-trained ion likelihood
    static constexpr double TRAIN_FDR = 0.05;

    /**
      @brief Training hits of self-trained models.

      The spectra are split into two halves by the parity of their "scan_index". Per half and class (cross-link or not,
      NuXL:isXL), the identifications whose best (first) hit is a target at q <= @p fdr, with q = (D + 1) / T over the
      hit score computed within the half (decoys win ties). Identifications without hits or scan index are skipped.

      @param[out] counts training hits per [half][class] (class 1: cross-links)
      @return indices into @p peptide_ids per half
    */
    static std::array<std::vector<Size>, 2> confidentHitsPerHalf(const PeptideIdentificationList& peptide_ids, double fdr,
                                                                 std::array<std::array<Size, 2>, 2>& counts);

    /// Weights of the native score (PXD048145, RNA-UV): see rescoreHits()
    static constexpr double DELTA_WEIGHT = 1.0;
    static constexpr double ION_LLR_WEIGHT = 0.3;

    /**
      @brief NuXL:rescoring native: score = NuXL:score + DELTA_WEIGHT * NuXL:delta_score (cross-links only)
      + ION_LLR_WEIGHT * NuXL:ion_llr, for every hit; then the hits of each identification are sorted and ranked again.

      The delta score is added to cross-links only: for non-cross-linked peptides it gained PSMs at 1% FDR but lost some
      at 10% FDR. Missing meta values count as 0.
    */
    static void rescoreHits(PeptideIdentificationList& peptide_ids);
  };
} // namespace OpenMS
