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

using namespace xercesc;
using namespace std;

namespace OpenMS
{

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
