// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathInferenceConfig.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathInferenceData.h>
#include <OpenMS/config.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Shared OpenSwath analysis context inference for peptide, protein, and gene rows.

    This class is file-format agnostic and operates only on compact typed rows.

    This is adapted from PyProphet's context inference methods. See the original publication for details: Rosenberger, G., Bludau, I., Schmitt, U. et al. Statistical control of peptide and protein error rates in large-scale targeted data-independent acquisition analyses. Nat Methods 14, 921–927 (2017). https://doi.org/10.1038/nmeth.4398

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI LevelContextInference
  {
  public:
    /**
      @brief Estimate p-values, q-values, and PEPs for the given input rows.

      The input rows are expected to contain one best-scoring row per entity and
      context. Their @p decoy flag defines the target/decoy split used to build
      score distributions for p-value, q-value, and local-FDR/PEP estimation.

      For @c run-specific inference, target/decoy error statistics are estimated
      independently per @c RUN_ID. For @c global and @c experiment-wide inference,
      one shared score distribution is used for all rows in the call.

      Internally this wraps OpenMS multiple-testing utilities for:
      - target/decoy-based p-value estimation
      - pi0 estimation
      - q-value estimation
      - local-FDR / posterior error probability estimation

      @param input Compact input rows selected for one inference level
      @param config Inference and error-estimation configuration

      @return Result rows containing the original entity and score plus derived
              p-value, q-value, and PEP estimates
    */
    static std::vector<LevelContextResultRow> infer(const std::vector<LevelContextInputRow>& input,
                                                    const LevelContextInferenceConfig& config);
  };
} // namespace OpenMS
