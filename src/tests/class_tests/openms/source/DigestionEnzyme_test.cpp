// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
///////////////////////////

#include <set>
#include <sstream>

using namespace OpenMS;
using namespace std;

// alias to force the regex constructor (the 2-arg form is ambiguous with the
// cut-specification constructor; a std::set 3rd argument selects this overload)
using SS = std::set<std::string>;

START_TEST(DigestionEnzyme, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DigestionEnzyme* ptr = nullptr;
DigestionEnzyme* null_ptr = nullptr;

START_SECTION((DigestionEnzyme(const std::string& name, const std::string& cleavage_regex, const std::set<std::string>& synonyms = std::set<std::string>(), std::string regex_description = "")))
{
  ptr = new DigestionEnzyme("Trypsin", "(?<=[KR])(?!P)", SS{"syn_a", "syn_b"}, "tryptic");
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_STRING_EQUAL(ptr->getName(), "Trypsin")
  TEST_STRING_EQUAL(ptr->getRegEx(), "(?<=[KR])(?!P)")
  TEST_STRING_EQUAL(ptr->getRegExDescription(), "tryptic")
  TEST_EQUAL(ptr->getSynonyms().size(), 2)
  TEST_EQUAL(ptr->getSynonyms().count("syn_a"), 1)
}
END_SECTION

START_SECTION((~DigestionEnzyme()))
{
  delete ptr;
}
END_SECTION

START_SECTION((DigestionEnzyme(const std::string& name, std::string cut_before, const std::string& nocut_after = "", std::string sense = "C", const std::set<std::string>& synonyms = std::set<std::string>(), std::string regex_description = "")))
{
  // the cut-specification constructor builds a cleavage regex
  DigestionEnzyme e("Tryp2", "KR", "P", "C");
  TEST_STRING_EQUAL(e.getName(), "Tryp2")
  TEST_NOT_EQUAL(e.getRegEx().empty(), true) // a regex was synthesized

  DigestionEnzyme e_n("AspN", "D", "", "N");
  TEST_NOT_EQUAL(e_n.getRegEx().empty(), true)

  // an empty cleavage specification is rejected
  TEST_EXCEPTION(Exception::MissingInformation, DigestionEnzyme("bad", "", "", "C"))
}
END_SECTION

START_SECTION((void setName(const std::string& name)))
{
  DigestionEnzyme e("E", "r", SS{});
  e.setName("Chymotrypsin");
  TEST_STRING_EQUAL(e.getName(), "Chymotrypsin")
}
END_SECTION

START_SECTION((const std::string& getName() const))
{
  DigestionEnzyme e("MyName", "r", SS{});
  TEST_STRING_EQUAL(e.getName(), "MyName")
}
END_SECTION

START_SECTION((void setSynonyms(const std::set<std::string>& synonyms)))
{
  DigestionEnzyme e("E", "r", SS{});
  e.setSynonyms(SS{"a", "b", "c"});
  TEST_EQUAL(e.getSynonyms().size(), 3)
  TEST_EQUAL(e.getSynonyms().count("b"), 1)
}
END_SECTION

START_SECTION((const std::set<std::string>& getSynonyms() const))
{
  DigestionEnzyme e("E", "r", SS{"only"});
  TEST_EQUAL(e.getSynonyms().size(), 1)
  TEST_EQUAL(e.getSynonyms().count("only"), 1)
}
END_SECTION

START_SECTION((void addSynonym(const std::string& synonym)))
{
  DigestionEnzyme e("E", "r", SS{});
  e.addSynonym("x");
  e.addSynonym("y");
  e.addSynonym("x"); // set semantics -> no duplicate
  TEST_EQUAL(e.getSynonyms().size(), 2)
}
END_SECTION

START_SECTION((void setRegEx(const std::string& cleavage_regex)))
{
  DigestionEnzyme e("E", "r", SS{});
  e.setRegEx("(?<=[FYW])");
  TEST_STRING_EQUAL(e.getRegEx(), "(?<=[FYW])")
}
END_SECTION

START_SECTION((const std::string& getRegEx() const))
{
  DigestionEnzyme e("E", "(?<=[KR])", SS{});
  TEST_STRING_EQUAL(e.getRegEx(), "(?<=[KR])")
}
END_SECTION

START_SECTION((void setRegExDescription(const std::string& value)))
{
  DigestionEnzyme e("E", "r", SS{});
  e.setRegExDescription("cuts after K/R");
  TEST_STRING_EQUAL(e.getRegExDescription(), "cuts after K/R")
}
END_SECTION

START_SECTION((const std::string& getRegExDescription() const))
{
  DigestionEnzyme e("E", "r", SS{}, "the description");
  TEST_STRING_EQUAL(e.getRegExDescription(), "the description")
}
END_SECTION

START_SECTION((bool operator==(const DigestionEnzyme& enzyme) const))
{
  DigestionEnzyme a("N", "R", SS{"s"}, "d");
  DigestionEnzyme b("N", "R", SS{"s"}, "d");
  DigestionEnzyme c("N2", "R", SS{"s"}, "d"); // different name
  TEST_EQUAL(a == b, true)
  TEST_EQUAL(a == c, false)
}
END_SECTION

START_SECTION((bool operator!=(const DigestionEnzyme& enzyme) const))
{
  DigestionEnzyme a("N", "R", SS{"s"}, "d");
  DigestionEnzyme c("N", "R2", SS{"s"}, "d"); // different regex
  TEST_EQUAL(a != c, true)
  TEST_EQUAL(a != a, false)
}
END_SECTION

START_SECTION((bool operator==(const std::string& cleavage_regex) const))
{
  DigestionEnzyme e("E", "(?<=[KR])(?!P)", SS{});
  TEST_EQUAL(e == std::string("(?<=[KR])(?!P)"), true)
  TEST_EQUAL(e == std::string("nope"), false)
}
END_SECTION

START_SECTION((bool operator!=(const std::string& cleavage_regex) const))
{
  DigestionEnzyme e("E", "(?<=[KR])(?!P)", SS{});
  TEST_EQUAL(e != std::string("other"), true)
  TEST_EQUAL(e != std::string("(?<=[KR])(?!P)"), false)
}
END_SECTION

START_SECTION((bool operator<(const DigestionEnzyme& enzyme) const))
{
  // ordering is by name
  DigestionEnzyme a("AAA", "r", SS{});
  DigestionEnzyme b("BBB", "r", SS{});
  TEST_EQUAL(a < b, true)
  TEST_EQUAL(b < a, false)
  TEST_EQUAL(a < a, false)
}
END_SECTION

START_SECTION((virtual bool setValueFromFile(const std::string& key, const std::string& value)))
{
  DigestionEnzyme e("E", "r", SS{});
  // recognized keys (matched by suffix / substring) return true and mutate
  TEST_EQUAL(e.setValueFromFile("Enzyme:Name", "NewName"), true)
  TEST_STRING_EQUAL(e.getName(), "NewName")

  TEST_EQUAL(e.setValueFromFile("Enzyme:RegEx", "(?<=[FY])"), true)
  TEST_STRING_EQUAL(e.getRegEx(), "(?<=[FY])")

  TEST_EQUAL(e.setValueFromFile("Enzyme:RegExDescription", "desc"), true)
  TEST_STRING_EQUAL(e.getRegExDescription(), "desc")

  TEST_EQUAL(e.setValueFromFile("Enzyme:Synonyms:1", "aka"), true)
  TEST_EQUAL(e.getSynonyms().count("aka"), 1)

  // an unrecognized key is ignored and returns false
  TEST_EQUAL(e.setValueFromFile("Enzyme:Bogus", "v"), false)
}
END_SECTION

START_SECTION([EXTRA] operator<<)
{
  DigestionEnzyme e("Trypsin", "(?<=[KR])(?!P)", SS{}, "tryptic");
  std::ostringstream os;
  os << e;
  TEST_STRING_EQUAL(os.str(), "digestion enzyme:Trypsin (cleavage: (?<=[KR])(?!P) - tryptic)")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
