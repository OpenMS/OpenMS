// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/UnimodXMLDataProvider.h>
#include <OpenMS/FORMAT/UnimodXMLFile.h>

using namespace std;

namespace OpenMS
{

  UnimodXMLDataProvider::UnimodXMLDataProvider(const String& filename)
    : filename_(filename)
  {
  }

  std::vector<std::unique_ptr<ResidueModification>> UnimodXMLDataProvider::loadModifications()
  {
    // UnimodXMLFile::load returns raw pointers; wrap them in unique_ptr
    vector<ResidueModification*> raw_mods;
    UnimodXMLFile().load(filename_, raw_mods);

    // Wrap all raw pointers first for exception safety, then mutate
    vector<unique_ptr<ResidueModification>> result;
    result.reserve(raw_mods.size());
    for (auto* m : raw_mods)
    {
      result.push_back(unique_ptr<ResidueModification>(m));
    }
    for (auto& m : result)
    {
      m->setFullId();
    }
    return result;
  }

} // namespace OpenMS
