// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Xiao Liang  $
// $Authors: Xiao Liang, Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/BuiltInProteaseDataProvider.h>
#include <OpenMS/CHEMISTRY/EnzymeXMLDataProvider.h>
#include <fstream>
using namespace std;

namespace OpenMS
{
  ProteaseDB::ProteaseDB():
    DigestionEnzymeDB<DigestionEnzymeProtein, ProteaseDB>()
  {
    // Create default providers: built-in enzymes + optional XML file for user overrides
    vector<unique_ptr<DigestionEnzymeDataProvider<DigestionEnzymeProtein>>> providers;
    providers.push_back(make_unique<BuiltInProteaseDataProvider>());
    providers.push_back(make_unique<EnzymeXMLDataProvider<DigestionEnzymeProtein>>("CHEMISTRY/Enzymes.xml", /*optional=*/true));
    loadFromProviders_(providers);
  }

  ProteaseDB::ProteaseDB(vector<unique_ptr<DigestionEnzymeDataProvider<DigestionEnzymeProtein>>> providers):
    DigestionEnzymeDB<DigestionEnzymeProtein, ProteaseDB>()
  {
    loadFromProviders_(providers);
  }

  void ProteaseDB::getAllXTandemNames(vector<String>& all_names) const
  {
    all_names.clear();
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if (!(*it)->getXTandemID().empty())
      {
        all_names.push_back((*it)->getName());
      }
    }
  }

  void ProteaseDB::getAllCometNames(vector<String>& all_names) const
  {
    all_names.clear();
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if ((*it)->getCometID() != -1)
      {
        all_names.push_back((*it)->getName());
      }
    }
  }

  void ProteaseDB::getAllOMSSANames(vector<String>& all_names) const
  {
    all_names.clear();
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if ((*it)->getOMSSAID() != -1)
      {
        all_names.push_back((*it)->getName());
      }
    }
  }

  void ProteaseDB::getAllMSGFNames(vector<String>& all_names) const
  {
    all_names.clear();
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if ((*it)->getMSGFID() != -1) // MS-GF+ starts enzyme numbering at 0
      {
        all_names.push_back((*it)->getName());
      }
    }
  }

  void ProteaseDB::writeTSV(String const& filename)
  {
    std::ofstream ofs(filename, std::ofstream::out);
    ofs << "OpenMS_AllowedEnzymes" << "\n";
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      ofs << (*it)->getName() << "\n";
    }
  }
} // namespace OpenMS
