// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/config.h>

#include <string>

namespace OpenMS
{
  /**
    @brief Helper Functions to perform an update query to the OpenMS REST server

    @ingroup System
  */
  class OPENMS_DLLAPI UpdateCheck
  {
public:
  static void run(const std::string& tool_name, const std::string& version, int debug_level);
  };
}


