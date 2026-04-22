// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Input to domain-agnostic Percolator::rescore.

    Row ordering is preserved in the output. Every row corresponds to one
    data point (PSM, transition, whatever the caller chooses).

    @var features         [n_rows][n_features] scalar features per row. Rows must
                          all have the same length.
    @var is_decoy         Target (false) or decoy (true) label per row.
    @var cv_group_keys    Per-row integer key used to group rows into the same
                          cross-validation fold. Rows sharing a key will never
                          be split across folds. Leave empty to use row index
                          (each row in its own group). Supply this when rows
                          have natural duplication (e.g., multiple PSMs from
                          one spectrum, multiple transitions from one precursor).
    @var feature_names    Names aligned 1:1 with feature columns; used for
                          logging only.
  */
  struct OPENMS_DLLAPI RescoreInput
  {
    std::vector<std::vector<double>> features;
    std::vector<bool> is_decoy;
    std::vector<int> cv_group_keys;
    StringList feature_names;

    /**
      @brief Optional per-row PIN-compatible fields.

      When supplied, these override the synthetic defaults
      (`scan = row_index`, `specFileNr = 0`, `expMass = calcMass = 0.0`) that
      the wrapper uses internally. Set them to make the in-process output
      bitwise-identical to running the external percolator binary on a .pin
      file derived from the same PSMs — Percolator's internal PSM sort order
      (OrderScanHash hashes `specFileNr` + `scan`) determines the CV fold
      assignment, so scan/spec_file_numbers have to match the PIN path.

      See Percolator::fillPINCompatibleFields() for a helper that populates
      these from a vector of PeptideIdentifications the way PercolatorInfile
      would write them.

      Each vector is either empty (use default) or has exactly n_rows entries
      (one per feature row).
    */
    std::vector<int>    scan_numbers;
    std::vector<int>    spec_file_numbers;
    std::vector<double> exp_masses;
    std::vector<double> calc_masses;
  };

  /// Output from Percolator::rescore. Aligned 1:1 with RescoreInput::features.
  struct OPENMS_DLLAPI RescoreOutput
  {
    std::vector<double> scores;     ///< SVM discriminant score per row
    std::vector<double> q_values;   ///< q-value per row
    std::vector<double> peps;       ///< posterior error probability per row
  };
}
