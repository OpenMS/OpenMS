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

using namespace xercesc;
using namespace std;

namespace OpenMS
{

  // Register UnimodXMLFile::load() as the Unimod loader callback so that
  // ModificationsDB (in Core) can load Unimod XML without depending on IO.
  static const bool unimod_loader_registered_ = []() {
    ModificationsDB::registerUnimodLoader([](const String& filename, std::vector<ResidueModification*>& mods) {
      UnimodXMLFile().load(filename, mods);
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
