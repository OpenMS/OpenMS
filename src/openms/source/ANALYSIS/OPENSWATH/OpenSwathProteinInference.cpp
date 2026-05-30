// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathProteinInference.h>

#include <OpenMS/ANALYSIS/OPENSWATH/LevelContextInference.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
  std::vector<LevelContextResultRow> OpenSwathProteinInference::infer(const std::vector<LevelContextInputRow>& input,
                                                                      const LevelContextInferenceConfig& config) const
  {
    if (config.level != InferenceLevel::Protein)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "OpenSwathProteinInference requires config.level = InferenceLevel::Protein.");
    }
    return LevelContextInference::infer(input, config);
  }
} // namespace OpenMS
