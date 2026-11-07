// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

namespace OpenMS::OpenDIA
{
  /** @brief A scored window with precursor and peptide-family provenance. */
  struct CNNQValueInput
  {
    std::string sample_id;
    std::string precursor_id;
    std::string group_id; ///< Peptide family, shared across related precursors.
    bool is_decoy = false;
    int fold = 0; ///< Identity of the held-out fold and its trained model.
    double score = 0.0; ///< Larger scores indicate stronger evidence.
    /// Optional target-decoy pair key for paired competition; empty means unpaired.
    std::string pair_key;
  };

  /** @brief Original window and its precursor-level selection and q-value. */
  struct CNNQValueResult
  {
    CNNQValueInput input;
    bool precursor_best = false;
    double q_value = 1.0; ///< Only precursor winners receive estimated q-values.
  };

  /**
    @brief Select precursor winners and estimate conservative q-values per fold.

    Select the highest score per precursor, breaking ties by lexical sample ID.
    Within each fold, group exact score ties and use cumulative target and decoy
    precursor counts T and D: min(1, (D+1)*ratio/T), or 1 when T is zero.
    Reverse cumulative minima yield q-values. Scores from different models are
    never pooled. All rows retain input order, with nonwinner q-values equal to 1.

    This is separate assay decoy-count estimation, not paired target-decoy
    competition. It guarantees neither FDR control nor prevalence correction.
    Empty input returns empty output, but the supplied ratio must still be valid.

    @param[in] inputs Scored windows with unique sample IDs, finite scores,
      nonnegative folds and nonblank IDs without tabs, carriage returns, newlines
      or embedded NUL bytes.
      Each precursor must have one label, group and fold, and each family one fold.
      Each represented fold must contain target and decoy winning precursors.
    @param[in] target_decoy_ratio Positive finite target-to-decoy assay ratio.
    @return One result per input window, in original order.
    @throws std::invalid_argument If input validation fails.
  */
  std::vector<CNNQValueResult> aggregateCNNQValues(
    const std::vector<CNNQValueInput>& inputs, double target_decoy_ratio = 1.0);

  /**
    @brief Paired (picked-precursor) target-decoy competition q-values per fold.

    Precursor winners are selected as in aggregateCNNQValues(). Winners sharing a
    nonempty pair_key form one target-decoy pair. Only the higher-scoring member
    competes (exact ties go to the decoy). Unpaired winners compete alone. Within
    each fold, competing entries use min(1, (D+1)/T) with exact-tie grouping and
    reverse cumulative minima. Pair members that lose keep q = 1.

    This is an experimental estimate, not guaranteed FDR control.

    @param[in] inputs As for aggregateCNNQValues(). A pair key may join at most one
      target and one decoy precursor, both in the same fold.
    @return One result per input window, in original order.
    @throws std::invalid_argument If input or pairing validation fails.
  */
  std::vector<CNNQValueResult> aggregatePairedCNNQValues(const std::vector<CNNQValueInput>& inputs);
}
