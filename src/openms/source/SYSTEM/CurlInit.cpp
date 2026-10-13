// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/CurlInit.h>

#ifdef WITH_CURL
#include <curl/curl.h>
#endif

namespace OpenMS
{
  void CurlInit::ensure()
  {
    static CurlInit instance;
    (void)instance;
  }

  // Without libcurl (WITH_CURL=OFF) there is nothing to initialize.
  CurlInit::CurlInit()
  {
#ifdef WITH_CURL
    curl_global_init(CURL_GLOBAL_DEFAULT);
#endif
  }

  CurlInit::~CurlInit()
  {
#ifdef WITH_CURL
    curl_global_cleanup();
#endif
  }
}
