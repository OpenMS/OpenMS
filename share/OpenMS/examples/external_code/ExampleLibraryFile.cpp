// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include "ExampleLibraryFile.h"

using namespace std;

namespace OpenMSExternal //optional namespace... however you like it
{
  std::string ExampleLibraryFile::printSomething()
  {
    return "this is the external library.";
  }

}
