// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/CurlInit.h>

#include <curl/curl.h>

namespace OpenMS
{
  void CurlInit::ensure()
  {
    static CurlInit instance;
    (void)instance;
  }

  CurlInit::CurlInit()
  {
    curl_global_init(CURL_GLOBAL_DEFAULT);
  }

  CurlInit::~CurlInit()
  {
    curl_global_cleanup();
  }
}
