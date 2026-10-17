// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Param.h>

namespace OpenMS
{
  /**
    @brief Shared transformation-model parameters for alignment algorithms and tools.

    Parameter values, constraints and descriptions follow the individual models.
  */
  class OPENMS_DLLAPI TransformationModelDefaults
  {
  public:
    /**
      @brief Assemble model selection and parameters for all supported models.

      @param[in] default_model Initial selection; an additional name (such as
      "none") is prepended to the allowed selections, preserving the tool convention.
    */
    static Param getDefaults(const std::string& default_model);
  };
} // namespace OpenMS
