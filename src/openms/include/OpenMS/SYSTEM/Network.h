// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <string>

namespace OpenMS
{
  /// Utility class for network operations (HTTP downloads)
  class OPENMS_DLLAPI Network
  {
  public:
    /**
      @brief Download file from given URL into a download folder. Returns when done.

      If a file with same filename already exists, creates a new file with '.number' appended to the basename.

      @throw IOException if download failed.
    */
    static void downloadFile(const std::string& url, const std::string& download_folder);
  };
}
