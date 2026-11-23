// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/CHEMISTRY/ModificationDataProvider.h>
#include <OpenMS/CHEMISTRY/UnimodXMLDataProvider.h>
#include <OpenMS/CHEMISTRY/OBODataProvider.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <iterator>

using namespace OpenMS;
using namespace std;

START_TEST(ModificationDataProvider, "$Id$")

/////////////////////////////////////////////////////////////
// InMemoryDataProvider tests
/////////////////////////////////////////////////////////////

START_SECTION(InMemoryDataProvider::loadModifications())
{
  vector<unique_ptr<ResidueModification>> test_mods;
  auto mod1 = make_unique<ResidueModification>();
  mod1->setId("TestMod1");
  mod1->setFullName("Test Modification 1");
  mod1->setOrigin('C');
  mod1->setDiffMonoMass(57.02);
  mod1->setFullId();
  test_mods.push_back(std::move(mod1));

  auto mod2 = make_unique<ResidueModification>();
  mod2->setId("TestMod2");
  mod2->setFullName("Test Modification 2");
  mod2->setOrigin('M');
  mod2->setDiffMonoMass(15.99);
  mod2->setFullId();
  test_mods.push_back(std::move(mod2));

  InMemoryDataProvider provider(std::move(test_mods));
  auto result = provider.loadModifications();

  TEST_EQUAL(result.size(), 2)
  TEST_STRING_EQUAL(result[0]->getId(), "TestMod1")
  TEST_STRING_EQUAL(result[1]->getId(), "TestMod2")
  TEST_REAL_SIMILAR(result[0]->getDiffMonoMass(), 57.02)
}
END_SECTION

/////////////////////////////////////////////////////////////
// UnimodXMLDataProvider tests
/////////////////////////////////////////////////////////////

START_SECTION(UnimodXMLDataProvider::loadModifications())
{
  UnimodXMLDataProvider provider("CHEMISTRY/unimod.xml");
  auto mods = provider.loadModifications();

  // unimod.xml should contain many modifications
  TEST_EQUAL(mods.size() > 100, true)

  // Check that FullId is set on all modifications
  bool all_have_fullid = true;
  for (const auto& m : mods)
  {
    if (m->getFullId().empty())
    {
      all_have_fullid = false;
      break;
    }
  }
  TEST_EQUAL(all_have_fullid, true)

  // Check that a well-known modification exists
  bool found_oxidation = false;
  for (const auto& m : mods)
  {
    if (m->getFullId() == "Oxidation (M)")
    {
      found_oxidation = true;
      break;
    }
  }
  TEST_EQUAL(found_oxidation, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// OBODataProvider tests (non-crosslink mode)
/////////////////////////////////////////////////////////////

START_SECTION(OBODataProvider::loadModifications() non-crosslink)
{
  OBODataProvider provider("CHEMISTRY/PSI-MOD.obo", false);
  auto mods = provider.loadModifications();

  // PSI-MOD should contain modifications
  TEST_EQUAL(mods.size() > 0, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// OBODataProvider tests (crosslink mode)
/////////////////////////////////////////////////////////////

START_SECTION(OBODataProvider::loadModifications() crosslink)
{
  OBODataProvider provider("CHEMISTRY/XLMOD.obo", true);
  auto mods = provider.loadModifications();

  // XLMOD in cross-link mode should contain cross-linkers
  TEST_EQUAL(mods.size() > 0, true)
}
END_SECTION

START_SECTION([EXTRA] OBODataProvider::loadModifications() reads a large file in parts as it reads it as a whole)
{
  // PSI-MOD.obo is read in parts that start at "[Term]" lines; with a blank after each "[Term]" (the same term
  // separator to the reader) no part can start, and the file is read as a whole
  std::string text;
  {
    std::ifstream in(File::find("CHEMISTRY/PSI-MOD.obo"), std::ios::binary);
    text.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
  }
  std::string whole_file;
  NEW_TMP_FILE(whole_file)
  {
    std::string whole;
    Size n_terms = 0;
    for (size_t from = 0; from < text.size();)
    {
      const size_t p = text.find("\n[Term]", from);
      if (p == std::string::npos) { whole += text.substr(from); break; }
      whole += text.substr(from, p + 7 - from) + " ";
      from = p + 7;
      ++n_terms;
    }
    TEST_EQUAL(n_terms > 1000, true)
    std::ofstream(whole_file, std::ios::binary) << whole;
  }
  for (bool cross_links_only : {false, true})
  {
    auto in_parts = OBODataProvider("CHEMISTRY/PSI-MOD.obo", cross_links_only).loadModifications();
    auto as_whole = OBODataProvider(whole_file, cross_links_only).loadModifications();
    TEST_EQUAL(in_parts.size(), as_whole.size())
    Size differences = 0;
    for (Size i = 0; i < std::min(in_parts.size(), as_whole.size()); ++i)
    {
      if (!(*in_parts[i] == *as_whole[i]) || in_parts[i]->getFullId() != as_whole[i]->getFullId() ||
          in_parts[i]->getSynonyms() != as_whole[i]->getSynonyms()) ++differences;
    }
    TEST_EQUAL(differences, 0)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// ModificationsDB with InMemoryDataProvider
/////////////////////////////////////////////////////////////

START_SECTION([EXTRA] ModificationsDB with InMemoryDataProvider)
{
  vector<unique_ptr<ResidueModification>> test_mods;
  auto mod = make_unique<ResidueModification>();
  mod->setId("Phospho");
  mod->setFullName("Phosphorylation");
  mod->setOrigin('S');
  mod->setDiffMonoMass(79.966);
  mod->setFullId();
  test_mods.push_back(std::move(mod));

  vector<unique_ptr<ModificationDataProvider>> providers;
  providers.push_back(make_unique<InMemoryDataProvider>(std::move(test_mods)));

  ModificationsDB db(std::move(providers));

  TEST_EQUAL(db.getNumberOfModifications(), 1)
  TEST_EQUAL(db.has("Phospho (S)"), true)
  TEST_EQUAL(db.has("Phospho"), true)
  TEST_EQUAL(db.has("Phosphorylation"), true)
  TEST_EQUAL(db.has("Nonexistent"), false)
  TEST_REAL_SIMILAR(db.getModification(0)->getDiffMonoMass(), 79.966)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
