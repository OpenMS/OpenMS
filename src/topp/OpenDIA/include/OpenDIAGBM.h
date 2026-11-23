// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenDIAScorer.h>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace OpenMS::OpenDIA
{
  /**
    @brief Deterministic histogram gradient boosting for binary logistic loss.

    Features are quantile binned per feature on the training rows (up to max_bins bins, uint8 codes). Trees grow
    leaf-wise up to max_leaf_nodes leaves; a split needs min_samples_leaf rows on each side. Split gain and leaf values
    use gradient and hessian sums with l2 regularisation (sums guarded by a small epsilon). The initial score is the
    weighted log-odds. Boosting stops early when the weighted validation logloss has not improved for `patience`
    rounds (all fitted trees are kept, as in sklearn). Histograms are built in parallel over features; every reduction
    has a fixed order, so the model does not depend on the thread count.
  */
  class HistGBM
  {
  public:
    struct Options
    {
      int max_iter = 300;
      double learning_rate = 0.05;
      int max_leaf_nodes = 31;
      int max_bins = 255;
      std::size_t min_samples_leaf = 20;
      double l2 = 0.0;
      int patience = 10;
      double tolerance = 1e-7;
      int threads = 1;
    };

    struct Node
    {
      int feature = -1;     ///< -1 for a leaf.
      float threshold = 0;  ///< Go left if x < threshold (NaN goes right).
      int left = -1, right = -1;
      double value = 0.0;   ///< Leaf output (shrinkage applied).
    };

    double init = 0.0;
    std::vector<std::vector<Node>> trees;

    /**
      @brief Fit on rows `train` of the row-major matrix `x` (n x d); rows `valid` drive early stopping (may be empty).
      @param label 1 positive, 0 negative, indexed like `x` rows.
      @param weight Sample weight, indexed like `x` rows.
    */
    void fit(const float* x, std::size_t d, const std::vector<std::uint32_t>& train, const std::vector<std::uint32_t>& valid,
             const std::vector<std::uint8_t>& label, const std::vector<double>& weight, const Options& options);

    /// Raw margin (log-odds) of one row.
    double predict(const float* row) const;
  };

  /**
    @brief Hard-negative GBM rescoring of CrossFitScorer windows.

    Up to `rounds` rounds (default 3): positives are the target winners at pooled paired TDC q <= `positive_q` (default 0.01) under the current score, negatives
    all rows of the top `hard_fraction` decoy precursors ranked by their winner score (0: all decoy winners). Per fold
    a HistGBM is trained on the other folds (negatives subsampled to at most 1.5 M, positive weight n_neg / n_pos,
    validation rows by pair_hash % 100 < 10) and scores its fold by the raw margin. Finally each precursor keeps its
    best window, scores are calibrated with CrossFitScorer::calibrateTail and q-values come from pooledPairedTDC.
    Only score and q_value of the result are filled.
  */
  CrossFitScorer::Result gbmRescore(const CrossFitScorer::Input& input, const std::vector<float>& score, double hard_fraction,
                                    std::uint64_t seed, int threads, const std::function<void(const std::string&)>& log = {},
                                    int rounds = 3, double positive_q = 0.01);
}
