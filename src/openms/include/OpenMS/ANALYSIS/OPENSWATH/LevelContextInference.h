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

    @ingroup OpenSwath
  */
  class OPENMS_DLLAPI LevelContextInference
  {
  public:
    /**
      @brief Estimate p-values, q-values, and PEPs for the given input rows.

      @param input Compact input rows selected for one inference level
      @param config Inference and error-estimation configuration
      @return Result rows with OpenSwath-style score lookups
    */
    static std::vector<LevelContextResultRow> infer(const std::vector<LevelContextInputRow>& input,
                                                    const LevelContextInferenceConfig& config);
  };
} // namespace OpenMS
