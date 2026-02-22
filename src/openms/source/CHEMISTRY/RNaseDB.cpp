// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser  $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/RNaseDB.h>

using namespace std;

namespace OpenMS
{
  RNaseDB::RNaseDB():
    DigestionEnzymeDB<DigestionEnzymeRNA, RNaseDB>()
  {
    vector<unique_ptr<DigestionEnzymeDataProvider<DigestionEnzymeRNA>>> providers;
    auto xml = createXmlProvider_("CHEMISTRY/Enzymes_RNA.xml", /*optional=*/false);
    if (xml) providers.push_back(std::move(xml));
    loadFromProviders_(providers);
  }

  RNaseDB::RNaseDB(vector<unique_ptr<DigestionEnzymeDataProvider<DigestionEnzymeRNA>>> providers):
    DigestionEnzymeDB<DigestionEnzymeRNA, RNaseDB>()
  {
    loadFromProviders_(providers);
  }

} // namespace OpenMS
