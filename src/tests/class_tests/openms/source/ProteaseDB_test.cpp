// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang, Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ProteaseDB, "$Id$")

/////////////////////////////////////////////////////////////

ProteaseDB* ptr = nullptr;
ProteaseDB* nullPointer = nullptr;
String RKP("(?<=[RX])(?!P)");
START_SECTION(ProteaseDB* getInstance())
    ptr = ProteaseDB::getInstance();
    TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(virtual ~ProteaseDB())
    NOT_TESTABLE
END_SECTION

START_SECTION((bool hasEnzyme(const String &name) const))
    TEST_EQUAL(ptr->hasEnzyme("Try"), false)
    TEST_EQUAL(ptr->hasEnzyme("Trypsin"), true)
END_SECTION

START_SECTION((const DigestionEnzymeProtein* getEnzyme(const String &name) const))
    TEST_EQUAL(ptr->getEnzyme("Trypsin")->getName(), "Trypsin")
    // test the synonyms
    TEST_EQUAL(ptr->getEnzyme("Clostripain")->getName(), "Arg-C")
    TEST_EXCEPTION(Exception::ElementNotFound, ptr->getEnzyme("DOESNOTEXIST"))
END_SECTION

START_SECTION((bool hasRegEx(const String& cleavage_regex) const))
    TEST_EQUAL(ptr->hasRegEx("(?<=[P])(?!P)"), false)
    TEST_EQUAL(ptr->hasRegEx(RKP), true)
END_SECTION

START_SECTION((const DigestionEnzymeProtein* getEnzymeByRegEx(const String& cleavage_regex) const))
    TEST_EQUAL(ptr->getEnzymeByRegEx(RKP)->getName(), "Arg-C")
END_SECTION

START_SECTION(bool hasEnzyme(const DigestionEnzymeProtein* enzyme) const)
  TEST_EQUAL(ptr->hasEnzyme(ptr->getEnzyme("Trypsin")), true)
  DigestionEnzymeProtein myNewEnzyme("bla", "blubb");
  TEST_EQUAL(ptr->hasEnzyme(&myNewEnzyme), false);
END_SECTION

START_SECTION(ConstEnzymeIterator beginEnzyme() const)
    ProteaseDB::EnzymeIterator it = ptr->beginEnzyme();
    Size count(0);
    while (it != ptr->endEnzyme())
    {
      ++it;
      ++count;
    }
    TEST_EQUAL(count >= 10, true)
END_SECTION

START_SECTION(ConstEnzymeIterator endEnzyme() const)
    NOT_TESTABLE // tested above
END_SECTION

START_SECTION((void getAllNames(std::vector<String>& all_names) const))
    vector<String> names;
    ptr->getAllNames(names);
    TEST_EQUAL(find(names.begin(), names.end(), "Trypsin") != names.end(), true)
    TEST_EQUAL(find(names.begin(), names.end(), "Tryptryp") != names.end(), false)
    Size old_size=names.size();
    ptr->getAllNames(names);
    TEST_EQUAL(names.size(), old_size)
END_SECTION

START_SECTION((void getAllXTandemNames(std::vector<String>& all_names) const))
    vector<String> names;
    ptr->getAllXTandemNames(names);
    TEST_EQUAL(find(names.begin(), names.end(), "Trypsin") != names.end(), true)
    TEST_EQUAL(find(names.begin(), names.end(), "no cleavage") != names.end(), true)
    Size old_size=names.size();
    ptr->getAllXTandemNames(names);
    TEST_EQUAL(names.size(), old_size)
END_SECTION

START_SECTION((void getAllOMSSANames(std::vector<String>& all_names) const))
    vector<String> names;
    ptr->getAllOMSSANames(names);
    TEST_EQUAL(find(names.begin(), names.end(), "Trypsin") != names.end(), true)
    TEST_EQUAL(find(names.begin(), names.end(), "leukocyte elastase") != names.end(), false)
    Size old_size=names.size();
    ptr->getAllOMSSANames(names);
    TEST_EQUAL(names.size(), old_size)
END_SECTION

START_SECTION([EXTRA] multithreaded example)
{

   int nr_iterations (1e2), test (0);
#pragma omp parallel for reduction (+: test)
  for (int k = 1; k < nr_iterations + 1; k++)
  {
    auto p = ProteaseDB::getInstance();
    int tmp (0);
    if (p->hasEnzyme("Trypsin"), true)
    {
      tmp++;
    }
    test += tmp;
  }
  TEST_EQUAL(test, nr_iterations)
}
END_SECTION

END_TEST
