// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>

#include <unordered_set>
#include <unordered_map>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(DigestionEnzyme, "$Id$")

/////////////////////////////////////////////////////////////

DigestionEnzyme* e_ptr = nullptr;
DigestionEnzyme* e_null = nullptr;

START_SECTION((DigestionEnzyme(const String& name, const String& cleavage_regex, const std::set<String>& synonyms, String regex_description)))
  set<String> syns;
  syns.insert("TestSyn1");
  syns.insert("TestSyn2");
  e_ptr = new DigestionEnzyme("TestEnzyme", "(?<=[RK])(?!P)", syns, "cuts after R or K but not before P");
  TEST_NOT_EQUAL(e_ptr, e_null)
  TEST_EQUAL(e_ptr->getName(), "TestEnzyme")
  TEST_EQUAL(e_ptr->getRegEx(), "(?<=[RK])(?!P)")
  TEST_EQUAL(e_ptr->getRegExDescription(), "cuts after R or K but not before P")
  TEST_EQUAL(e_ptr->getSynonyms().size(), 2)
END_SECTION

START_SECTION((virtual ~DigestionEnzyme()))
  delete e_ptr;
END_SECTION

set<String> syns;
syns.insert("Syn1");
e_ptr = new DigestionEnzyme("Enzyme1", "(?<=[RK])", syns, "cuts after R or K");

START_SECTION(DigestionEnzyme(const DigestionEnzyme& enzyme))
  DigestionEnzyme copy(*e_ptr);
  TEST_EQUAL(copy, *e_ptr)
END_SECTION

START_SECTION(DigestionEnzyme& operator=(const DigestionEnzyme& enzyme))
  DigestionEnzyme copy("", "", std::set<String>());
  copy = *e_ptr;
  TEST_EQUAL(copy, *e_ptr)
END_SECTION

START_SECTION(void setName(const String& name))
  DigestionEnzyme copy(*e_ptr);
  e_ptr->setName("NewName");
  TEST_NOT_EQUAL(copy, *e_ptr)
END_SECTION

START_SECTION(const String& getName() const)
  TEST_EQUAL(e_ptr->getName(), "NewName")
END_SECTION

START_SECTION(void setSynonyms(const std::set<String>& synonyms))
  DigestionEnzyme copy(*e_ptr);
  set<String> new_syns;
  new_syns.insert("NewSyn");
  e_ptr->setSynonyms(new_syns);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(void addSynonym(const String& synonym))
  DigestionEnzyme copy(*e_ptr);
  e_ptr->addSynonym("AnotherSyn");
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(const std::set<String>& getSynonyms() const)
  TEST_EQUAL(e_ptr->getSynonyms().size(), 2)
END_SECTION

START_SECTION(void setRegEx(const String& cleavage_regex))
  DigestionEnzyme copy(*e_ptr);
  e_ptr->setRegEx("(?<=[P])");
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(const String& getRegEx() const)
  TEST_EQUAL(e_ptr->getRegEx(), "(?<=[P])")
END_SECTION

START_SECTION(void setRegExDescription(const String& value))
  DigestionEnzyme copy(*e_ptr);
  e_ptr->setRegExDescription("new description");
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(const String& getRegExDescription() const)
  TEST_EQUAL(e_ptr->getRegExDescription(), "new description")
END_SECTION

START_SECTION(bool operator==(const DigestionEnzyme& enzyme) const)
  DigestionEnzyme r("", "", std::set<String>());
  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  r.setName("other_name");
  TEST_EQUAL(r == *e_ptr, false)

  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  r.setRegEx("?<=[X]");
  TEST_EQUAL(r == *e_ptr, false)

  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  set<String> new_syns;
  new_syns.insert("different_syn");
  r.setSynonyms(new_syns);
  TEST_EQUAL(r == *e_ptr, false)

  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  r.setRegExDescription("different description");
  TEST_EQUAL(r == *e_ptr, false)
END_SECTION

START_SECTION(bool operator!=(const DigestionEnzyme& enzyme) const)
  DigestionEnzyme r("", "", std::set<String>());
  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setName("other_name");
  TEST_EQUAL(r != *e_ptr, true)
END_SECTION

START_SECTION(bool operator==(const String& cleavage_regex) const)
  TEST_EQUAL(*e_ptr == "(?<=[P])", true)
  TEST_EQUAL(*e_ptr == "(?<=[X])", false)
END_SECTION

START_SECTION(bool operator!=(const String& cleavage_regex) const)
  TEST_EQUAL(*e_ptr != "(?<=[X])", true)
  TEST_EQUAL(*e_ptr != "(?<=[P])", false)
END_SECTION

START_SECTION(bool operator<(const DigestionEnzyme& enzyme) const)
  DigestionEnzyme a("AAA", "", std::set<String>());
  DigestionEnzyme b("BBB", "", std::set<String>());
  TEST_EQUAL(a < b, true)
  TEST_EQUAL(b < a, false)
END_SECTION

START_SECTION([EXTRA] std::hash<DigestionEnzyme>)
{
  // Test that equal objects have equal hashes
  set<String> syns1;
  syns1.insert("Syn1");
  syns1.insert("Syn2");
  DigestionEnzyme e1("TestEnzyme", "(?<=[RK])(?!P)", syns1, "test description");
  DigestionEnzyme e2("TestEnzyme", "(?<=[RK])(?!P)", syns1, "test description");

  std::hash<DigestionEnzyme> hasher;
  TEST_EQUAL(hasher(e1), hasher(e2))

  // Test that different objects (likely) have different hashes
  DigestionEnzyme e3("DifferentEnzyme", "(?<=[RK])(?!P)", syns1, "test description");
  TEST_NOT_EQUAL(hasher(e1), hasher(e3))

  // Test use in unordered_set
  std::unordered_set<DigestionEnzyme> enzyme_set;
  enzyme_set.insert(e1);
  enzyme_set.insert(e2); // Should not increase size (duplicate)
  enzyme_set.insert(e3);
  TEST_EQUAL(enzyme_set.size(), 2)
  TEST_EQUAL(enzyme_set.count(e1), 1)
  TEST_EQUAL(enzyme_set.count(e3), 1)

  // Test use in unordered_map
  std::unordered_map<DigestionEnzyme, int> enzyme_map;
  enzyme_map[e1] = 1;
  enzyme_map[e2] = 2; // Should overwrite e1's value
  enzyme_map[e3] = 3;
  TEST_EQUAL(enzyme_map.size(), 2)
  TEST_EQUAL(enzyme_map[e1], 2)
  TEST_EQUAL(enzyme_map[e3], 3)

  // Test that hash changes when fields change
  DigestionEnzyme e4 = e1;
  std::size_t hash_before = hasher(e4);
  e4.setName("ChangedName");
  std::size_t hash_after = hasher(e4);
  TEST_NOT_EQUAL(hash_before, hash_after)

  // Test with different synonyms
  set<String> syns2;
  syns2.insert("DifferentSyn");
  DigestionEnzyme e5("TestEnzyme", "(?<=[RK])(?!P)", syns2, "test description");
  TEST_NOT_EQUAL(hasher(e1), hasher(e5))
}
END_SECTION

delete e_ptr;

END_TEST
