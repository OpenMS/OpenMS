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
  /// Thread-safe RAII guard for curl_global_init / curl_global_cleanup.
  /// Call CurlInit::ensure() before using any curl API.
  class CurlInit
  {
  public:
    OPENMS_DLLAPI static void ensure();

  private:
    CurlInit();
    ~CurlInit();
    CurlInit(const CurlInit&) = delete;
    CurlInit& operator=(const CurlInit&) = delete;
  };
}
