// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include "ExampleLibraryFile.h"

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

using namespace std;
using namespace OpenMS;

//optional namespace... however you like it
namespace OpenMSExternal
{
  std::string ExampleLibraryFile::printSomething()
  {
    return "this is the external library.";
  }

  void ExampleLibraryFile::loadAndSaveFeatureXML()
  {
    FeatureMap fm;
    Feature feature;
    fm.push_back(feature);
    String tmpfilename = "tmpfile.featureXML";
    FeatureXMLFile().store(tmpfilename, fm);

    FeatureMap fm2;
    FeatureXMLFile().store(tmpfilename, fm2);
  }
}
