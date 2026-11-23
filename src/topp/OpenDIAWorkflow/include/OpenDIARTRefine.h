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
#include <string>
#include <vector>

namespace OpenMS::OpenDIA
{
  /**
    @brief Per-run RT refinement (default since PLAN_RTR): a ridge regression of the RT residual on sequence composition.

    Features per precursor: token counts over a vocabulary of modified residues, N- and C-terminal residue one-hots,
    residue count and its log, a charge one-hot (1..4) and the precursor m/z. Features are standardised on the
    training rows; the ridge penalty is chosen by 5-fold cross-validation with folds given by the caller.
  */
  class RTRefineModel
  {
  public:
    /// Tokens of a modified sequence: a residue letter with an optional "(UniMod:N)" directly after it, or a bare
    /// "(UniMod:N)" (N-terminal modification). Leading '.' characters are skipped.
    static std::vector<std::string> tokenize(const std::string& sequence);

    /// Build the vocabulary (sorted distinct tokens) and the residue alphabet from all library sequences.
    explicit RTRefineModel(const std::vector<std::string>& sequences);

    /// Feature vector length.
    std::size_t dimension() const { return vocabulary_.size() + 2 * alphabet_.size() + 2 + 4 + 1; }

    /// Features of one precursor (tokens outside the vocabulary are ignored).
    std::vector<double> features(const std::string& sequence, int charge, double mz) const;

    const std::vector<std::string>& vocabulary() const { return vocabulary_; }
    const std::vector<char>& alphabet() const { return alphabet_; }

  private:
    std::vector<std::string> vocabulary_;
    std::vector<char> alphabet_;
  };

  /// Solve A x = b for a symmetric positive definite d x d matrix A (row major) by Cholesky; false if not SPD.
  bool choleskySolve(std::vector<double> a, std::vector<double> b, std::size_t d, std::vector<double>& x);

  /// A fitted ridge model with intercept on standardised features.
  struct RidgeFit
  {
    std::vector<double> mean, sd, weights;
    double intercept = 0.0;
    double predict(const double* row) const;
  };

  /// Fit ridge on the rows @p rows of the n x d matrix @p x (row major); @p y[k] is the target of row rows[k].
  RidgeFit fitRidge(const std::vector<double>& x, std::size_t d, const std::vector<std::size_t>& rows, const std::vector<double>& y, double alpha);

  struct RTRefineResult
  {
    std::vector<double> correction;     ///< Final-model prediction for every row of x (seconds).
    std::vector<double> out_of_fold;    ///< Cross-fitted predictions for the training rows, in training order.
    double alpha = 0.0;
    std::size_t training = 0;
  };

  /**
    @brief Choose alpha from 10^{-2, -1.5, ..., 4} by 5-fold CV (minimum mean squared error, first on ties), then fit.

    @param[in] x         n x d features (row major), all precursors
    @param[in] d         Feature count
    @param[in] train     Training row indices into x
    @param[in] y         Targets, one per training row
    @param[in] folds     Fold (0..4) per training row
  */
  RTRefineResult fitRTRefine(const std::vector<double>& x, std::size_t d, const std::vector<std::size_t>& train, const std::vector<double>& y,
                             const std::vector<std::uint32_t>& folds);
}
