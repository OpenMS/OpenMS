// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/TypeAliases.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Input to domain-agnostic Percolator::rescore.

    Row ordering is preserved in the output. Each row corresponds to one
    data point — the row's semantics (PSM, transition, peak group, etc.)
    are determined by the caller.
  */
  struct OPENMS_DLLAPI RescoreInput
  {
    /// [n_rows][n_features] scalar features per row. Rows must all have
    /// the same length. Empty if the rows are in the dense @ref feature_matrix.
    std::vector<std::vector<double>> features;
    /// Target (false) or decoy (true) label per row.
    std::vector<bool> is_decoy;
    /// Per-row integer key used to group rows into the same cross-validation
    /// fold. Rows sharing a key will never be split across folds. Leave empty
    /// to use row index (each row in its own group). Supply this when rows
    /// have natural duplication (e.g., multiple PSMs from one spectrum,
    /// multiple transitions from one precursor).
    std::vector<int> cv_group_keys;
    /// Names aligned 1:1 with feature columns; used for logging only.
    StringList feature_names;

    /**
      @brief Optional per-row PIN-compatible fields.

      When supplied, these override the synthetic defaults
      (`scan = row_index`, `specFileNr = 0`, `expMass = calcMass = 0.0`)
      that the wrapper uses internally. Populate them when the in-process
      output must be bitwise-identical to running the external percolator
      binary on a .pin file derived from the same PSMs: Percolator's
      internal PSM sort order (`OrderScanHash` hashes `specFileNr` and
      `scan`) determines the CV fold assignment, so scan_numbers and
      spec_file_numbers must match those that `PercolatorInfile::store`
      would emit.

      Percolator::fillPINCompatibleFields() is a helper that derives all
      four vectors from a vector of PeptideIdentifications using the same
      conventions as PercolatorInfile::store.

      Each vector is either empty (in which case the default is used) or
      contains exactly n_rows entries, one per feature row.
    */
    std::vector<int>    scan_numbers;
    std::vector<int>    spec_file_numbers;
    std::vector<double> exp_masses;
    std::vector<double> calc_masses;

    /**
      @brief The feature matrix in one dense block, in place of @ref features.

      Row-major: row i holds the values of the feature_names.size() features at
      [i * feature_names.size(), (i + 1) * feature_names.size()). Used when
      @ref features is empty (at most one of the two may hold rows); the number
      of rows is then is_decoy.size(). rowCount(), featureCount() and row()
      read the rows of either layout.

      Percolator::buildPSMRescoreInput() builds the feature matrix of PSM
      rescoring this way, straight from the PSMs.
    */
    std::vector<double> feature_matrix;

    /// Number of rows: of @ref features, or else of the dense @ref feature_matrix
    size_t rowCount() const
    {
      return features.empty() && !feature_matrix.empty() ? is_decoy.size() : features.size();
    }

    /// Number of features per row: of the first row of @ref features, or else feature_names.size() (0 without rows)
    size_t featureCount() const
    {
      if (!features.empty()) return features.front().size();
      return feature_matrix.empty() ? 0 : feature_names.size();
    }

    /// The featureCount() values of row @p i: of @ref features, or else of the dense @ref feature_matrix
    const double* row(size_t i) const
    {
      return !features.empty() ? features[i].data() : feature_matrix.data() + i * feature_names.size();
    }
  };

  /// Output from Percolator::rescore. Aligned 1:1 with the rows of the RescoreInput (RescoreInput::rowCount()).
  struct OPENMS_DLLAPI RescoreOutput
  {
    std::vector<double> scores;     ///< SVM discriminant score per row
    std::vector<double> q_values;   ///< q-value per row
    std::vector<double> peps;       ///< posterior error probability per row
  };

  /**
    @brief Trained Percolator model: averaged SVM weights in raw feature space.

    Produced by Percolator::train and consumed by Percolator::score. The
    weights are un-normalized: they are intended to multiply raw input
    features directly. The normalization transform learned by the SVM has
    already been folded into the weights and bias by
    Normalizer::unnormalizeweight(). Callers must therefore not normalize
    features before score(); doing so would apply the transform twice.

    The raw SVM dot product for a row with feature vector f is
        raw = sum_j(f[j] * weights[j]) + weights[n_features]     // bias last
    Percolator::score() applies a further FDR-based rescaling on top of
    this raw value to produce the final SVM discriminant reported in
    RescoreOutput.scores; see Percolator::score for the exact formula.
  */
  struct OPENMS_DLLAPI PercolatorModel
  {
    /// Integer schema version for the on-disk format.
    int format_version = 1;
    /// Feature column names. Must be non-empty and must match
    /// RescoreInput::feature_names positionally at score time. Any
    /// string value is permitted: the bias is stored in the header by
    /// saveModel(), so feature names carry no reserved meaning.
    StringList feature_names;
    /// Size = n_features + 1. The last entry is the bias.
    std::vector<double> weights;
    /// "stdv" | "uni" | "none" — the normalizer used during training.
    /// Informational only; all three produce raw-space weights that
    /// score() can apply directly, since the normalization transform is
    /// already folded into the weights and bias. Recorded so that
    /// reproducibility tooling can identify the learner configuration
    /// that produced the model.
    std::string normalizer_type;
    /// Random seed used during training. Informational.
    int seed = 0;
  };
}
