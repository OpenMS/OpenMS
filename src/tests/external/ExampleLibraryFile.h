// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <string>

//optional namespace... however you like it
namespace OpenMSExternal
{

  class ExampleLibraryFile
  {
public:
    static std::string printSomething();

    // just to have a dependency to OpenMS in the lib
    void loadAndSaveFeatureXML();
  };

}

