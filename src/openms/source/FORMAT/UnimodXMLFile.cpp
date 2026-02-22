// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/UnimodXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/UnimodXMLHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/UnimodXMLDataProvider.h>

using namespace xercesc;
using namespace std;

namespace OpenMS
{

  // Register a UnimodXMLDataProvider factory so that ModificationsDB (in Core)
  // can load Unimod XML without depending on IO types at link time.
  static const bool unimod_provider_registered_ = []() {
    ModificationsDB::registerXmlProviderFactory([](const String& filename) {
      return std::make_unique<UnimodXMLDataProvider>(filename);
    });
    return true;
  }();

  UnimodXMLFile::UnimodXMLFile() :
    Internal::XMLFile()
  {

  }

  UnimodXMLFile::~UnimodXMLFile() = default;

  void UnimodXMLFile::load(const String& filename, vector<ResidueModification*> & modifications)
  {
    String file = File::find(filename);

    Internal::UnimodXMLHandler handler(modifications, file);
    parse_(file, &handler);
  }

} // namespace OpenMS
