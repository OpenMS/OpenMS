// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace OpenMS::OpenDIA
{
  /**
    @brief Small network whose output is scored by a decoy-whitened linear readout.

    For a standardised feature vector z (dimension d) the embedding is h = [z, tanh(W2 relu(W1 z + b1) + b2)] with
    `hidden` and `embed` units. For positive rows P and negative rows N the readout is
    w = (S_N + lambda I)^-1 (mu_P - mu_N), with S_N the sample covariance of the negatives, and its objective
    obj = (mu_P - mu_N)' w is the squared Mahalanobis separation in decoy units. Training maximises obj.
    Weights and biases start uniform in +-1/sqrt(fan_in).
  */
  class WhitenedNetwork
  {
  public:
    WhitenedNetwork(std::size_t dimension, std::size_t hidden, std::size_t embed, std::uint64_t seed);

    std::size_t dimension() const { return dimension_; }
    /// Embedding size d + embed.
    std::size_t outputs() const { return dimension_ + embed_; }
    std::vector<double>& parameters() { return parameters_; }
    const std::vector<double>& parameters() const { return parameters_; }

    /// h of one row (outputs() values).
    void embed(const float* z, double* h) const;

    /**
      @brief Readout objective of a batch and, optionally, its gradient with respect to parameters().
      @param[in] positive Feature rows of the positives (at least 1).
      @param[in] negative Feature rows of the negatives (at least 2).
      @param[in] lambda Ridge added to the negative covariance (> 0).
      @param[out] gradient If not null, d obj / d parameters (resized).
      @param[in] threads Worker threads; the result does not depend on it.
    */
    double objective(const std::vector<const float*>& positive, const std::vector<const float*>& negative, double lambda,
                     std::vector<double>* gradient, int threads) const;

    /**
      @brief Readout from embeddings: w = (S_N + lambda I)^-1 (mu_P - mu_N).
      @param[in] positive Row-major positive embeddings (k values per row).
      @param[in] negative Row-major negative embeddings (at least 2 rows).
      @param[in] k Embedding size.
      @param[in] lambda Ridge (> 0).
      @param[out] objective (mu_P - mu_N)' w.
      @return w.
    */
    static std::vector<double> readout(const std::vector<double>& positive, const std::vector<double>& negative, std::size_t k,
                                       double lambda, double& objective);

  private:
    std::size_t dimension_, hidden_, embed_;
    std::vector<double> parameters_; ///< W1 [hidden x d], b1, W2 [embed x hidden], b2.
  };

  /**
    @brief Cross-fitted, self-trained scorer for candidate windows with pooled paired target-decoy q-values.

    Rows are candidate windows (peaks) of precursors. Each precursor and its target-decoy pair belong to one fold. For
    each fold k a model is trained on the other folds only and scores fold k:

    1. Features are standardised by the mean and SD of the training folds.
    2. Round 0 trains on the seed rows (e.g. anchor windows as positives, decoy windows as negatives).
    3. Rounds 1 and 2 relabel: the training folds are scored by the previous model, precursor winners compete in the
       pooled paired TDC, target winners at q <= self_q[r] are the positives and all decoy precursor winners the negatives.
    4. Each fit uses Adam on minibatches of batch_positive + batch_negative rows for at most max_epochs epochs of
       steps_per_epoch steps, and keeps the parameters with the best readout objective on a validation share of the
       rows (selected by a hash of the row's group, so pairs or families stay together; patience `patience` epochs).
       The final readout is computed on the non-validation rows, and the score is h' w.

    Then, per precursor, only the window with the largest linear part z' w_z of its fold's model is kept (the W rule;
    others get NOT_KEPT). Scores are calibrated per fold to decoy tail units: s' = -log S_k(s), with S_k the decoy winner
    survival (interpolated, (count + 1) / (n + 1)) up to the top 1 % decoy score t_k, and an exponential tail above it
    whose scale is the mean excess of those decoys. Last, the pooled paired TDC gives q-values.

    No labels other than target / decoy are used. Results depend on the seed but not on the thread count.
  */
  class CrossFitScorer
  {
  public:
    static constexpr float NOT_KEPT = -1e30f;

    struct Options
    {
      std::uint64_t seed = 0;
      int threads = 1;
      std::size_t hidden = 64;
      std::size_t embed = 16;
      double lambda = 1e-2;
      std::size_t batch_positive = 256;
      std::size_t batch_negative = 1024;
      int max_epochs = 30;
      int steps_per_epoch = 100;
      int patience = 5;
      double learning_rate = 1e-3;
      unsigned validation_percent = 15; ///< Rows whose group hash % 100 is below this validate.
      std::vector<double> self_q = {0.15, 0.05}; ///< One self-training round per entry.
    };

    struct Seed
    {
      std::uint32_t row = 0;
      bool positive = false;
      std::uint64_t group_hash = 0; ///< Validation split key (e.g. of the peptide family).
    };

    struct Input
    {
      std::size_t dimension = 0;
      std::vector<float> features;            ///< rows x dimension, row-major.
      std::vector<std::uint32_t> precursor;   ///< Dense precursor index per row.
      std::vector<std::uint32_t> pair;        ///< Pair index per row; a target and its decoy precursor share it.
      std::vector<std::uint8_t> decoy;        ///< 1 for decoy rows.
      std::vector<std::uint8_t> fold;         ///< Fold per row; one per pair.
      std::vector<std::uint64_t> pair_hash;   ///< Validation split key per row (same within a pair).
      std::vector<Seed> seeds;                ///< Round-0 training rows (any fold; each fold's model uses the others).
      std::size_t rows() const { return precursor.size(); }
    };

    struct Result
    {
      std::vector<float> score;     ///< Calibrated score of each precursor's kept window, NOT_KEPT otherwise.
      std::vector<float> raw_score; ///< h' w of the row's fold model.
      std::vector<float> linear;    ///< z' w_z of the row's fold model.
      std::vector<double> q_value;  ///< Pooled paired TDC q of competing precursor winners, 1 otherwise.
      std::vector<std::string> log;
    };

    /**
      @brief Train, score and compute q-values.
      @throws std::invalid_argument On inconsistent input (sizes, folds within a pair, a fold without decoys).
    */
    static Result run(const Input& input, const Options& options, const std::function<void(const std::string&)>& log = {});

    /// Pooled paired TDC result over a subset of rows.
    struct Competition
    {
      std::vector<std::uint32_t> winners;   ///< Best row of each precursor present in the subset.
      std::vector<std::uint32_t> ranked;    ///< The better member of each pair, by standardised score, descending.
      std::vector<double> q;                ///< q per entry of `ranked`.
    };

    /**
      @brief Pooled paired target-decoy competition.

      Precursor winners are the highest-scoring row of each precursor (ties: lower row). Scores of winners are
      standardised per fold by the median and 1.4826 MAD of that fold's decoy winners. The better winner of each pair
      competes (ties: the decoy). Over all folds together q = min(1, (D + 1) / T) with a reverse cumulative minimum.
      @throws std::invalid_argument If a fold has no decoy winners or a zero decoy spread.
    */
    static Competition pooledPairedTDC(const Input& input, const std::vector<float>& score, const std::vector<std::uint32_t>& rows);

    /**
      @brief Per-fold decoy-tail calibration of the kept windows (rows at NOT_KEPT stay there).
      @throws std::invalid_argument If a fold has fewer than 10 decoy winners.
    */
    static std::vector<float> calibrateTail(const Input& input, const std::vector<float>& score);

    /// FNV-1a 64 of a string (the validation split hash).
    static std::uint64_t fnv1a(const std::string& s);

    /**
      @brief W rule, tail calibration and pooled paired TDC over the allowed rows only (score, q_value of @p res
      recomputed from raw_score and linear; rows not allowed are NOT_KEPT).
    */
    static void finalise(const Input& in, Result& res, const std::vector<bool>& allowed);

    /**
      @brief Choose the RT window half width per fold from the run itself (cross-fitted) and apply it.

      For fold k, every grid width w = fraction x @p tolerance is tried on the other folds: rows with
      |rt_error| > w are dropped, finalise() is rerun, and target precursors at q <= 0.01 are counted. The width with
      the most targets is fold k's (ties: the larger width). Then every fold is restricted to its own width and
      finalise() gives the reported scores and q-values. Uses only target/decoy status; scores are not retrained.
      @param[in] rt_error Per row, apex minus (calibrated) library RT in seconds.
      @param[in] fractions Grid of fractions of @p tolerance, ascending.
      @return The width of each fold.
    */
    static std::vector<double> selectWidth(const Input& in, Result& res, const std::vector<double>& rt_error, double tolerance,
                                           const std::vector<double>& fractions);
  };
}
