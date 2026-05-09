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
    @brief OpenSwath protein-level context inference.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI OpenSwathProteinInference
  {
  public:
    /// Infer protein-level statistics for the provided rows.
    std::vector<LevelContextResultRow> infer(const std::vector<LevelContextInputRow>& input,
                                             const LevelContextInferenceConfig& config) const;
  };
} // namespace OpenMS
