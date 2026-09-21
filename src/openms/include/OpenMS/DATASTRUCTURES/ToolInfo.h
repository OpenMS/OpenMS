// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Ruben Grünberg, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief Tool metadata shared by parameter-file serializers.
  */
  struct ToolInfo
  {
    std::string version_;
    std::string name_;
    std::string docurl_;
    std::string category_;
    std::string description_;
    std::vector<std::string> citations_;
  };
}
