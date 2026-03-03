// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/CHEMISTRY/RibonucleotideDataProvider.h>
#include <OpenMS/CHEMISTRY/ModomicsJSONDataProvider.h>
#include <OpenMS/CHEMISTRY/RibonucleotideTSVDataProvider.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CHEMISTRY/Ribonucleotide.h>

using namespace OpenMS;
using namespace std;

START_TEST(RibonucleotideDataProvider, "$Id$")

/////////////////////////////////////////////////////////////
// InMemoryRibonucleotideDataProvider tests
/////////////////////////////////////////////////////////////

START_SECTION(InMemoryRibonucleotideDataProvider::loadRibonucleotides())
{
  vector<RibonucleotideEntry> entries;

  RibonucleotideEntry e1;
  e1.ribo = make_unique<Ribonucleotide>();
  e1.ribo->setCode("TestA");
  e1.ribo->setName("Test Adenosine");
  entries.push_back(std::move(e1));

  RibonucleotideEntry e2;
  e2.ribo = make_unique<Ribonucleotide>();
  e2.ribo->setCode("TestC");
  e2.ribo->setName("Test Cytidine");
  entries.push_back(std::move(e2));

  InMemoryRibonucleotideDataProvider provider(std::move(entries));
  auto result = provider.loadRibonucleotides();

  TEST_EQUAL(result.size(), 2)
  TEST_STRING_EQUAL(result[0].ribo->getCode(), "TestA")
  TEST_STRING_EQUAL(result[0].ribo->getName(), "Test Adenosine")
  TEST_STRING_EQUAL(result[1].ribo->getCode(), "TestC")
  TEST_STRING_EQUAL(result[1].ribo->getName(), "Test Cytidine")
  TEST_EQUAL(result[0].isAmbiguous(), false)
  TEST_EQUAL(result[1].isAmbiguous(), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ModomicsJSONDataProvider tests
/////////////////////////////////////////////////////////////

START_SECTION(ModomicsJSONDataProvider::loadRibonucleotides())
{
  ModomicsJSONDataProvider provider("CHEMISTRY/Modomics.json");
  auto entries = provider.loadRibonucleotides();

  // Modomics JSON should contain many ribonucleotides
  TEST_EQUAL(entries.size() > 100, true)

  // Check that a well-known entry exists: Am = 2'-O-methyladenosine
  bool found_am = false;
  for (const auto& e : entries)
  {
    if (e.ribo->getCode() == "Am")
    {
      found_am = true;
      TEST_STRING_EQUAL(e.ribo->getName(), "2'-O-methyladenosine")
      break;
    }
  }
  TEST_EQUAL(found_am, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// RibonucleotideTSVDataProvider tests
/////////////////////////////////////////////////////////////

START_SECTION(RibonucleotideTSVDataProvider::loadRibonucleotides())
{
  RibonucleotideTSVDataProvider provider("CHEMISTRY/Custom_RNA_modifications.tsv");
  auto entries = provider.loadRibonucleotides();

  // TSV file should contain entries
  TEST_EQUAL(entries.size() > 0, true)

  // Check that the ambiguous entry "msU?" exists
  bool found_msu = false;
  for (const auto& e : entries)
  {
    if (e.ribo->getCode() == "msU?")
    {
      found_msu = true;
      TEST_EQUAL(e.isAmbiguous(), true)
      break;
    }
  }
  TEST_EQUAL(found_msu, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// RibonucleotideDB with InMemoryRibonucleotideDataProvider
/////////////////////////////////////////////////////////////

START_SECTION([EXTRA] RibonucleotideDB with InMemoryRibonucleotideDataProvider)
{
  vector<RibonucleotideEntry> entries;

  RibonucleotideEntry e1;
  e1.ribo = make_unique<Ribonucleotide>();
  e1.ribo->setCode("xA");
  e1.ribo->setName("Test Modified Adenosine");
  entries.push_back(std::move(e1));

  RibonucleotideEntry e2;
  e2.ribo = make_unique<Ribonucleotide>();
  e2.ribo->setCode("xC");
  e2.ribo->setName("Test Modified Cytidine");
  entries.push_back(std::move(e2));

  vector<unique_ptr<RibonucleotideDataProvider>> providers;
  providers.push_back(make_unique<InMemoryRibonucleotideDataProvider>(std::move(entries)));

  RibonucleotideDB db(std::move(providers));

  // verify getRibonucleotide works
  auto ribo_a = db.getRibonucleotide("xA");
  TEST_STRING_EQUAL(ribo_a->getCode(), "xA")
  TEST_STRING_EQUAL(ribo_a->getName(), "Test Modified Adenosine")

  auto ribo_c = db.getRibonucleotide("xC");
  TEST_STRING_EQUAL(ribo_c->getCode(), "xC")
  TEST_STRING_EQUAL(ribo_c->getName(), "Test Modified Cytidine")

  // verify iterator count
  Size count = 0;
  for (auto it = db.begin(); it != db.end(); ++it)
  {
    ++count;
  }
  TEST_EQUAL(count, 2)

  // verify ElementNotFound exception for unknown code
  TEST_EXCEPTION(Exception::ElementNotFound, db.getRibonucleotide("Nonexistent"))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
