// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/METADATA/ProteinModificationSummary.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

///////////////////////////

START_TEST(ProteinModificationSummary, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ProteinModificationSummary* ptr = nullptr;
ProteinModificationSummary* null_ptr = nullptr;

START_SECTION((ProteinModificationSummary()))
{
  ptr = new ProteinModificationSummary();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->AALevelSummary.empty(), true)
}
END_SECTION

START_SECTION((~ProteinModificationSummary()))
{
  delete ptr;
}
END_SECTION

START_SECTION(([ProteinModificationSummary::Statistics] bool operator==(const Statistics& rhs) const))
{
  ProteinModificationSummary::Statistics a;
  // default-constructed values
  TEST_EQUAL(a.count, 0)
  TEST_REAL_SIMILAR(a.frequency, 0.0)
  TEST_REAL_SIMILAR(a.FLR, 0.0)
  TEST_REAL_SIMILAR(a.probability, 0.0)

  ProteinModificationSummary::Statistics b;
  TEST_EQUAL(a == b, true)

  a.count = 123;
  a.frequency = 0.5;
  a.FLR = 0.01;
  a.probability = 0.99;
  TEST_EQUAL(a == b, false)

  b.count = 123;
  b.frequency = 0.5;
  b.FLR = 0.01;
  b.probability = 0.99;
  TEST_EQUAL(a == b, true)

  // each field participates in the comparison
  b.count = 124;
  TEST_EQUAL(a == b, false)
  b.count = 123;
  b.probability = 0.98;
  TEST_EQUAL(a == b, false)
}
END_SECTION

START_SECTION((bool operator==(const ProteinModificationSummary& rhs) const))
{
  const ResidueModification* ox = ModificationsDB::getInstance()->getModification("Oxidation", "M", ResidueModification::ANYWHERE);
  TEST_NOT_EQUAL(ox, 0)

  ProteinModificationSummary::Statistics st;
  st.count = 42;
  st.frequency = 0.25;
  st.FLR = 0.05;
  st.probability = 0.9;

  ProteinModificationSummary s1, s2;
  TEST_EQUAL(s1 == s2, true)

  // position 10 carries an Oxidation observed in 42 PSMs
  s1.AALevelSummary[10][ox] = st;
  TEST_EQUAL(s1 == s2, false)

  s2.AALevelSummary[10][ox] = st;
  TEST_EQUAL(s1 == s2, true)

  // differing statistic at the same position/modification
  ProteinModificationSummary::Statistics st2 = st;
  st2.count = 43;
  s2.AALevelSummary[10][ox] = st2;
  TEST_EQUAL(s1 == s2, false)

  // differing position
  s2.AALevelSummary[10][ox] = st;
  TEST_EQUAL(s1 == s2, true)
  s2.AALevelSummary[11][ox] = st;
  TEST_EQUAL(s1 == s2, false)
}
END_SECTION

START_SECTION([EXTRA] (nested type aliases and map semantics))
{
  const ResidueModification* ox = ModificationsDB::getInstance()->getModification("Oxidation", "M", ResidueModification::ANYWHERE);

  ProteinModificationSummary s;
  ProteinModificationSummary::Statistics st;
  st.count = 7;
  s.AALevelSummary[3][ox] = st;

  // ModificationsToStatistics: modification -> statistic
  ProteinModificationSummary::ModificationsToStatistics& at3 = s.AALevelSummary[3];
  TEST_EQUAL(at3.size(), 1)
  TEST_EQUAL(at3[ox].count, 7)

  // AALevelSummary: position -> ModificationsToStatistics
  TEST_EQUAL(s.AALevelSummary.size(), 1)
  TEST_EQUAL(s.AALevelSummary.count(3), 1)
  TEST_EQUAL(s.AALevelSummary.count(99), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
