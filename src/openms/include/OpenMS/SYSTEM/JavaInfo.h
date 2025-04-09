// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

namespace OpenMS
{
  class String;
  /**
    @brief Detect Java and retrieve information.

    Similar classes exist for other external tools, e.g. PythonInfo .

    @ingroup System
  */
  class OPENMS_DLLAPI JavaInfo
  {
public:
    /**
      @brief Determine if Java is installed and reachable

      The call fails if either Java is not installed or if a relative location is given and Java is not on the search PATH.

      @param java_executable Path to Java executable. Can be absolute, relative or just a filename
      @param verbose_on_error On error, should an error message be printed to OPENMS_LOG_ERROR?
      @return Returns false if Java executable can not be called; true if Java executable can be executed
    **/
    static bool canRun(const String& java_executable, bool verbose_on_error = true);
  };

}

