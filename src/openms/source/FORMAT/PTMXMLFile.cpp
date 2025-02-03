// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PTMXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/PTMXMLHandler.h>

using namespace std;

namespace OpenMS
{

  PTMXMLFile::PTMXMLFile() = default;

  void PTMXMLFile::load(const String & filename, map<String, pair<String, String> > & ptm_informations)
  {
    ptm_informations.clear();

    Internal::PTMXMLHandler handler(ptm_informations, filename);
    parse_(filename, &handler);
  }

  void PTMXMLFile::store(const String& filename, map<String, pair<String, String> > & ptm_informations) const
  {
    Internal::PTMXMLHandler handler(ptm_informations, filename);
    save_(filename, &handler);
  }

} // namespace OpenMS
