// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/config.h>

namespace OpenMS
{
  class String;
  /**
    @brief Helper Functions to perform an update query to the OpenMS REST server

    @ingroup System
  */
  class OPENMS_DLLAPI UpdateCheck
  {
public:
  static void run(const String& tool_name, const String& version, int debug_level);
  };
}


