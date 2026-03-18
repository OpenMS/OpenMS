// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/CrossLinksDB.h>
#include <OpenMS/CHEMISTRY/OBODataProvider.h>

using namespace std;

namespace OpenMS
{
  CrossLinksDB::CrossLinksDB()
    : ModificationsDB(makeCrossLinkProviders_())
  {
  }

  std::vector<std::unique_ptr<ModificationDataProvider>> CrossLinksDB::makeCrossLinkProviders_()
  {
    std::vector<std::unique_ptr<ModificationDataProvider>> providers;
    providers.push_back(std::make_unique<OBODataProvider>("CHEMISTRY/XLMOD.obo", /*cross_links_only=*/true));
    return providers;
  }


  CrossLinksDB::~CrossLinksDB() = default;

  void CrossLinksDB::getAllSearchModifications(vector<String>& modifications) const
  {
    modifications.clear();

    for (vector<ResidueModification*>::const_iterator it = mods_.begin(); it != mods_.end(); ++it)
    {
      if (!(*it)->getPSIMODAccession().empty())
      {
        modifications.push_back((*it)->getFullId());
      }
    }
    sort(modifications.begin(), modifications.end());
  }

}
