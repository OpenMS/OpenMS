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
#include <OpenMS/SYSTEM/File.h>

#include <filesystem>
#include <stdexcept>

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
#ifdef OPENMS_EXPECTED_DATA_DIR
    if (!std::filesystem::equivalent(File::getOpenMSDataPath(), OPENMS_EXPECTED_DATA_DIR))
    {
      throw std::runtime_error("OpenMS did not use the installed runtime data directory");
    }
#endif
    FeatureMap fm;
    Feature feature;
    fm.push_back(feature);
    std::string tmpfilename = "tmpfile.featureXML";
    FileHandler().storeFeatures(tmpfilename, fm, {FileTypes::FEATUREXML});

    FeatureMap fm2;
    FileHandler().loadFeatures(tmpfilename, fm2, {FileTypes::FEATUREXML});
    if (fm2.size() != fm.size())
    {
      throw std::runtime_error("Installed consumer featureXML round trip lost features");
    }
  }
}
