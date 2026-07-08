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
#include <OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DigestionEnzymeRNA, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DigestionEnzymeRNA* ptr = nullptr;
DigestionEnzymeRNA* null_ptr = nullptr;

START_SECTION([EXTRA] DigestionEnzymeRNA() and RNase typedef)
{
  // the base default constructor is inherited implicitly; RNA-specific fields
  // start out empty
  ptr = new DigestionEnzymeRNA();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_STRING_EQUAL(ptr->getCutsAfterRegEx(), "")
  TEST_STRING_EQUAL(ptr->getThreePrimeGain(), "")
  // RNase is an alias for DigestionEnzymeRNA
  RNase r;
  r.setName("RNase_T1");
  TEST_STRING_EQUAL(r.getName(), "RNase_T1")
  delete ptr;
}
END_SECTION

START_SECTION((void setCutsAfterRegEx(const std::string& value)))
{
  DigestionEnzymeRNA e;
  e.setCutsAfterRegEx("G");
  TEST_STRING_EQUAL(e.getCutsAfterRegEx(), "G")
}
END_SECTION

START_SECTION((std::string getCutsAfterRegEx() const))
{
  DigestionEnzymeRNA e;
  TEST_STRING_EQUAL(e.getCutsAfterRegEx(), "")
  e.setCutsAfterRegEx("AC");
  TEST_STRING_EQUAL(e.getCutsAfterRegEx(), "AC")
}
END_SECTION

START_SECTION((void setCutsBeforeRegEx(const std::string& value)))
{
  DigestionEnzymeRNA e;
  e.setCutsBeforeRegEx("U");
  TEST_STRING_EQUAL(e.getCutsBeforeRegEx(), "U")
}
END_SECTION

START_SECTION((std::string getCutsBeforeRegEx() const))
{
  DigestionEnzymeRNA e;
  TEST_STRING_EQUAL(e.getCutsBeforeRegEx(), "")
  e.setCutsBeforeRegEx("GA");
  TEST_STRING_EQUAL(e.getCutsBeforeRegEx(), "GA")
}
END_SECTION

START_SECTION((void setThreePrimeGain(const std::string& value)))
{
  DigestionEnzymeRNA e;
  e.setThreePrimeGain("p"); // 3' phosphate
  TEST_STRING_EQUAL(e.getThreePrimeGain(), "p")
}
END_SECTION

START_SECTION((std::string getThreePrimeGain() const))
{
  DigestionEnzymeRNA e;
  TEST_STRING_EQUAL(e.getThreePrimeGain(), "")
  e.setThreePrimeGain("OH");
  TEST_STRING_EQUAL(e.getThreePrimeGain(), "OH")
}
END_SECTION

START_SECTION((void setFivePrimeGain(const std::string& value)))
{
  DigestionEnzymeRNA e;
  e.setFivePrimeGain("OH"); // 5' hydroxyl
  TEST_STRING_EQUAL(e.getFivePrimeGain(), "OH")
}
END_SECTION

START_SECTION((std::string getFivePrimeGain() const))
{
  DigestionEnzymeRNA e;
  TEST_STRING_EQUAL(e.getFivePrimeGain(), "")
  e.setFivePrimeGain("p");
  TEST_STRING_EQUAL(e.getFivePrimeGain(), "p")
}
END_SECTION

START_SECTION((bool setValueFromFile(const std::string& key, const std::string& value) override))
{
  DigestionEnzymeRNA e;

  // RNA-specific keys (matched by suffix) set the corresponding member
  TEST_EQUAL(e.setValueFromFile("RNase:CutsAfter", "G"), true)
  TEST_STRING_EQUAL(e.getCutsAfterRegEx(), "G")
  TEST_EQUAL(e.setValueFromFile("RNase:CutsBefore", "A"), true)
  TEST_STRING_EQUAL(e.getCutsBeforeRegEx(), "A")
  TEST_EQUAL(e.setValueFromFile("RNase:ThreePrimeGain", "p"), true)
  TEST_STRING_EQUAL(e.getThreePrimeGain(), "p")
  TEST_EQUAL(e.setValueFromFile("RNase:FivePrimeGain", "OH"), true)
  TEST_STRING_EQUAL(e.getFivePrimeGain(), "OH")

  // base-class keys are still handled (delegated to DigestionEnzyme)
  TEST_EQUAL(e.setValueFromFile("RNase:Name", "RNase_X"), true)
  TEST_STRING_EQUAL(e.getName(), "RNase_X")

  // an unrecognized key returns false
  TEST_EQUAL(e.setValueFromFile("RNase:Bogus", "v"), false)
}
END_SECTION

START_SECTION([EXTRA] inherited operator== compares only the base-class fields)
{
  // operator== is inherited from DigestionEnzyme and compares name/regex/
  // synonyms/description only -- the RNA-specific fields are NOT part of it.
  DigestionEnzymeRNA a;
  a.setName("RNase"); a.setRegEx("(?<=G)"); a.setThreePrimeGain("p");
  DigestionEnzymeRNA b;
  b.setName("RNase"); b.setRegEx("(?<=G)"); b.setThreePrimeGain("different");
  TEST_EQUAL(a == b, true)  // equal despite differing 3' gains

  DigestionEnzymeRNA c;
  c.setName("Other"); c.setRegEx("(?<=G)");
  TEST_EQUAL(a == c, false) // differ in the base name
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
