// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Init.h>

#include <xercesc/util/PlatformUtils.hpp>

#include <mutex>

namespace OpenMS::Internal
{
  namespace
  {
    std::mutex& xercesMutex()
    {
      static std::mutex mutex;
      return mutex;
    }
  }

  void xercesInitialize()
  {
    const std::lock_guard<std::mutex> lock(xercesMutex());
    xercesc::XMLPlatformUtils::Initialize();
  }

  void xercesTerminate()
  {
    const std::lock_guard<std::mutex> lock(xercesMutex());
    xercesc::XMLPlatformUtils::Terminate();
  }

  // Initialize xerces
  // see ticket #352 for more details
  struct xerces_init
  {
    xerces_init() 
    {
      xercesInitialize();
    }

    ~xerces_init() 
    {
      xercesTerminate();
    }

  };
  const xerces_init xinit;

} //OpenMS //Internal

