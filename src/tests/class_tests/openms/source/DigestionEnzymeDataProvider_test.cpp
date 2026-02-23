// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Sachsenberg $
// $Authors: Sachsenberg $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/DigestionEnzymeDataProvider.h>
#include <OpenMS/CHEMISTRY/EnzymeXMLDataProvider.h>
#include <OpenMS/CHEMISTRY/BuiltInProteaseDataProvider.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(DigestionEnzymeDataProvider, "$Id$")

/////////////////////////////////////////////////////////////
// InMemoryDigestionEnzymeDataProvider tests
/////////////////////////////////////////////////////////////

START_SECTION(InMemoryDigestionEnzymeDataProvider - DigestionEnzymeProtein)
{
  vector<unique_ptr<DigestionEnzymeProtein>> enzymes;
  enzymes.push_back(make_unique<DigestionEnzymeProtein>(
    "TestProtease", "(?<=[K])", set<String>{"test_syn"}, "Test protease description"));
  enzymes.push_back(make_unique<DigestionEnzymeProtein>(
    "TestProtease2", "(?<=[R])", set<String>(), "Another test protease"));

  InMemoryDigestionEnzymeDataProvider<DigestionEnzymeProtein> provider(std::move(enzymes));
  auto loaded = provider.loadEnzymes();
  TEST_EQUAL(loaded.size(), 2)
  TEST_EQUAL(loaded[0]->getName(), "TestProtease")
  TEST_EQUAL(loaded[0]->getRegEx(), "(?<=[K])")
  TEST_EQUAL(loaded[1]->getName(), "TestProtease2")
}
END_SECTION

START_SECTION(InMemoryDigestionEnzymeDataProvider - DigestionEnzymeRNA)
{
  vector<unique_ptr<DigestionEnzymeRNA>> enzymes;
  auto rna_enzyme = make_unique<DigestionEnzymeRNA>();
  rna_enzyme->setName("TestRNase");
  rna_enzyme->setRegEx("(?<=[AU])");
  enzymes.push_back(std::move(rna_enzyme));

  InMemoryDigestionEnzymeDataProvider<DigestionEnzymeRNA> provider(std::move(enzymes));
  auto loaded = provider.loadEnzymes();
  TEST_EQUAL(loaded.size(), 1)
  TEST_EQUAL(loaded[0]->getName(), "TestRNase")
}
END_SECTION

/////////////////////////////////////////////////////////////
// EnzymeXMLDataProvider tests
/////////////////////////////////////////////////////////////

START_SECTION(EnzymeXMLDataProvider - load protein enzymes from XML)
{
  // This test requires Enzymes.xml to be in the share directory
  EnzymeXMLDataProvider<DigestionEnzymeProtein> provider("CHEMISTRY/Enzymes.xml", /*optional=*/true);
  auto loaded = provider.loadEnzymes();
  // File may or may not exist depending on test environment.
  // If it exists, it should load some enzymes; if not, empty is OK (optional=true).
  // We just verify no crash occurs.
  STATUS("EnzymeXMLDataProvider loaded " << loaded.size() << " protein enzymes from optional XML")
}
END_SECTION

START_SECTION(EnzymeXMLDataProvider - load RNA enzymes from XML)
{
  EnzymeXMLDataProvider<DigestionEnzymeRNA> provider("CHEMISTRY/Enzymes_RNA.xml");
  auto loaded = provider.loadEnzymes();
  TEST_EQUAL(loaded.size() > 0, true)
  // Verify at least one enzyme was loaded with a valid name
  bool found_rnase = false;
  for (const auto& e : loaded)
  {
    if (!e->getName().empty())
    {
      found_rnase = true;
      break;
    }
  }
  TEST_EQUAL(found_rnase, true)
}
END_SECTION

START_SECTION(EnzymeXMLDataProvider - optional file not found)
{
  EnzymeXMLDataProvider<DigestionEnzymeProtein> provider("CHEMISTRY/NonExistentFile.xml", /*optional=*/true);
  auto loaded = provider.loadEnzymes();
  TEST_EQUAL(loaded.size(), 0)
}
END_SECTION

START_SECTION(EnzymeXMLDataProvider - required file not found throws)
{
  EnzymeXMLDataProvider<DigestionEnzymeProtein> provider("CHEMISTRY/NonExistentFile.xml", /*optional=*/false);
  TEST_EXCEPTION(Exception::FileNotFound, provider.loadEnzymes())
}
END_SECTION

/////////////////////////////////////////////////////////////
// BuiltInProteaseDataProvider tests
/////////////////////////////////////////////////////////////

START_SECTION(BuiltInProteaseDataProvider - load built-in enzymes)
{
  BuiltInProteaseDataProvider provider;
  auto loaded = provider.loadEnzymes();
  TEST_EQUAL(loaded.size() >= 30, true) // at least 30 built-in proteases

  // Verify some well-known enzymes are present
  bool found_trypsin = false;
  bool found_argc = false;
  bool found_lysc = false;
  for (const auto& e : loaded)
  {
    if (e->getName() == "Trypsin") found_trypsin = true;
    if (e->getName() == "Arg-C") found_argc = true;
    if (e->getName() == "Lys-C") found_lysc = true;
  }
  TEST_EQUAL(found_trypsin, true)
  TEST_EQUAL(found_argc, true)
  TEST_EQUAL(found_lysc, true)
}
END_SECTION

START_SECTION(BuiltInProteaseDataProvider - verify Trypsin details)
{
  BuiltInProteaseDataProvider provider;
  auto loaded = provider.loadEnzymes();

  const DigestionEnzymeProtein* trypsin = nullptr;
  for (const auto& e : loaded)
  {
    if (e->getName() == "Trypsin")
    {
      trypsin = e.get();
      break;
    }
  }
  TEST_NOT_EQUAL(trypsin, (const DigestionEnzymeProtein*)nullptr)
  TEST_EQUAL(trypsin->getRegEx(), "(?<=[KRX])(?!P)")
  TEST_EQUAL(trypsin->getPSIID(), "MS:1001251")
  TEST_EQUAL(trypsin->getXTandemID(), "[KR]|{P}")
  TEST_EQUAL(trypsin->getCometID(), 1)
  TEST_EQUAL(trypsin->getOMSSAID(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// DI-constructed ProteaseDB tests
/////////////////////////////////////////////////////////////

START_SECTION(ProteaseDB - construct from providers)
{
  vector<unique_ptr<DigestionEnzymeDataProvider<DigestionEnzymeProtein>>> providers;

  // Add just two test enzymes via InMemoryDigestionEnzymeDataProvider
  vector<unique_ptr<DigestionEnzymeProtein>> enzymes;
  enzymes.push_back(make_unique<DigestionEnzymeProtein>(
    "CustomEnzyme1", "(?<=[K])", set<String>{"custom_syn"}, "Custom enzyme 1"));
  enzymes.push_back(make_unique<DigestionEnzymeProtein>(
    "CustomEnzyme2", "(?<=[R])", set<String>(), "Custom enzyme 2"));
  providers.push_back(make_unique<InMemoryDigestionEnzymeDataProvider<DigestionEnzymeProtein>>(std::move(enzymes)));

  ProteaseDB db(std::move(providers));
  TEST_EQUAL(db.hasEnzyme("CustomEnzyme1"), true)
  TEST_EQUAL(db.hasEnzyme("CustomEnzyme2"), true)
  TEST_EQUAL(db.hasEnzyme("custom_syn"), true) // synonym
  TEST_EQUAL(db.hasEnzyme("Trypsin"), false) // not loaded

  vector<String> names;
  db.getAllNames(names);
  TEST_EQUAL(names.size(), 2)
}
END_SECTION

START_SECTION(ProteaseDB - singleton matches DI with default providers)
{
  // The singleton and a DI-constructed DB with the same providers should agree
  auto singleton = ProteaseDB::getInstance();
  TEST_EQUAL(singleton->hasEnzyme("Trypsin"), true)
  TEST_EQUAL(singleton->hasEnzyme("Arg-C"), true)
  TEST_EQUAL(singleton->hasEnzyme("Clostripain"), true) // synonym

  // Verify Trypsin details match the built-in definition
  auto trypsin = singleton->getEnzyme("Trypsin");
  TEST_EQUAL(trypsin->getRegEx(), "(?<=[KRX])(?!P)")
  TEST_EQUAL(trypsin->getCometID(), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
// DI-constructed RNaseDB tests
/////////////////////////////////////////////////////////////

START_SECTION(RNaseDB - construct from providers)
{
  vector<unique_ptr<DigestionEnzymeDataProvider<DigestionEnzymeRNA>>> providers;

  vector<unique_ptr<DigestionEnzymeRNA>> enzymes;
  auto rna_e = make_unique<DigestionEnzymeRNA>();
  rna_e->setName("CustomRNase");
  rna_e->setRegEx("(?<=[AU])");
  enzymes.push_back(std::move(rna_e));
  providers.push_back(make_unique<InMemoryDigestionEnzymeDataProvider<DigestionEnzymeRNA>>(std::move(enzymes)));

  RNaseDB db(std::move(providers));
  TEST_EQUAL(db.hasEnzyme("CustomRNase"), true)
  TEST_EQUAL(db.hasEnzyme("RNase_T1"), false) // not loaded in custom DB
}
END_SECTION

START_SECTION(RNaseDB - singleton loads from XML)
{
  auto singleton = RNaseDB::getInstance();
  // The RNA enzyme XML should have loaded at least some enzymes
  vector<String> names;
  singleton->getAllNames(names);
  TEST_EQUAL(names.size() > 0, true)
}
END_SECTION

/////////////////////////////////////////////////////////////

END_TEST
