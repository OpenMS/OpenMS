// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser  $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/RNaseDB.h>
#include <OpenMS/CHEMISTRY/EnzymeXMLDataProvider.h>

using namespace std;

namespace OpenMS
{
  RNaseDB::RNaseDB():
    DigestionEnzymeDB<DigestionEnzymeRNA, RNaseDB>()
  {
    vector<unique_ptr<DigestionEnzymeDataProvider<DigestionEnzymeRNA>>> providers;
    providers.push_back(make_unique<EnzymeXMLDataProvider<DigestionEnzymeRNA>>("CHEMISTRY/Enzymes_RNA.xml"));
    loadFromProviders_(providers);
  }

  RNaseDB::RNaseDB(vector<unique_ptr<DigestionEnzymeDataProvider<DigestionEnzymeRNA>>> providers):
    DigestionEnzymeDB<DigestionEnzymeRNA, RNaseDB>()
  {
    loadFromProviders_(providers);
  }

} // namespace OpenMS
