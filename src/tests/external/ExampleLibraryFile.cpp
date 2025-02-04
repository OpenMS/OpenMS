// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include "ExampleLibraryFile.h"

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FileHandler.h>

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
    FileHandler().storeFeatures(tmpfilename, fm, {FileTypes::FEATUREXML});

    FeatureMap fm2;
    FileHandler().storeFeatures(tmpfilename, fm2, {FileTypes::FEATUREXML});
  }
}
