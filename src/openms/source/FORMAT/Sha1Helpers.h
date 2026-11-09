// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <string>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief SHA-1 (RFC 3174) of the contents of the file @p filename, as 40 lowercase hex characters.

      Private helper shared by FileHandler::computeFileHash and the imzML writer. Not installed.

      @return The hex digest, or an empty string if the file cannot be opened.
    */
    std::string sha1FileHex(const std::string& filename);
  }
}
