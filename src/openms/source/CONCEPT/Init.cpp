// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Init.h>

#include <xercesc/util/PlatformUtils.hpp>

namespace OpenMS::Internal
{

  // Initialize xerces
  // see ticket #352 for more details
  struct xerces_init
  {
    xerces_init() 
    {
      xercesc::XMLPlatformUtils::Initialize();
    }

    ~xerces_init() 
    {
      xercesc::XMLPlatformUtils::Terminate();
    }

  };
  const xerces_init xinit;

} //OpenMS //Internal

