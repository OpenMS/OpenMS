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

#include <OpenMS/ANALYSIS/NUXL/NuXLFragmentAdductDefinition.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <unordered_set>

///////////////////////////

START_TEST(NuXLFragmentAdductDefinition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

NuXLFragmentAdductDefinition* ptr = nullptr;
NuXLFragmentAdductDefinition* null_ptr = nullptr;

START_SECTION((NuXLFragmentAdductDefinition()))
{
  ptr = new NuXLFragmentAdductDefinition();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_STRING_EQUAL(ptr->name, "")
  TEST_REAL_SIMILAR(ptr->mass, 0.0)
  TEST_EQUAL(ptr->formula.isEmpty(), true)
}
END_SECTION

START_SECTION((~NuXLFragmentAdductDefinition()))
{
  delete ptr;
}
END_SECTION

START_SECTION((NuXLFragmentAdductDefinition(const EmpiricalFormula& f, const std::string& n, double m)))
{
  NuXLFragmentAdductDefinition fad(EmpiricalFormula("C5H10"), "adductA", 70.078);
  TEST_STRING_EQUAL(fad.name, "adductA")
  TEST_REAL_SIMILAR(fad.mass, 70.078)
  TEST_EQUAL(fad.formula == EmpiricalFormula("C5H10"), true)
}
END_SECTION

START_SECTION((NuXLFragmentAdductDefinition(const NuXLFragmentAdductDefinition&)))
{
  NuXLFragmentAdductDefinition fad(EmpiricalFormula("C5H10"), "adductA", 70.078);
  NuXLFragmentAdductDefinition copy(fad);
  TEST_STRING_EQUAL(copy.name, "adductA")
  TEST_REAL_SIMILAR(copy.mass, 70.078)
  TEST_EQUAL(copy.formula == EmpiricalFormula("C5H10"), true)
}
END_SECTION

START_SECTION((bool operator==(const NuXLFragmentAdductDefinition& other) const))
{
  EmpiricalFormula f("C5H10");
  NuXLFragmentAdductDefinition a(f, "adductA", 70.0);

  // equality is defined over (formula, name) only -- the mass is intentionally ignored
  NuXLFragmentAdductDefinition same_diff_mass(f, "adductA", 999.0);
  TEST_EQUAL(a == same_diff_mass, true)

  // different formula or name -> not equal
  NuXLFragmentAdductDefinition diff_formula(EmpiricalFormula("C6H12"), "adductA", 70.0);
  TEST_EQUAL(a == diff_formula, false)
  NuXLFragmentAdductDefinition diff_name(f, "adductD", 70.0);
  TEST_EQUAL(a == diff_name, false)

  // identical definition compares equal
  NuXLFragmentAdductDefinition copy(f, "adductA", 70.0);
  TEST_EQUAL(a == copy, true)
}
END_SECTION

START_SECTION((bool operator<(const NuXLFragmentAdductDefinition& other) const))
{
  EmpiricalFormula f("C5H10");

  // ordering is primarily by mass
  NuXLFragmentAdductDefinition lo(f, "x", 10.0);
  NuXLFragmentAdductDefinition hi(f, "x", 20.0);
  TEST_EQUAL(lo < hi, true)
  TEST_EQUAL(hi < lo, false)

  // note: two definitions equal under operator== (same formula+name) can still be
  // ordered by operator< if their masses differ (operator< also considers mass)
  NuXLFragmentAdductDefinition a(f, "adductA", 70.0);
  NuXLFragmentAdductDefinition b(f, "adductA", 999.0);
  TEST_EQUAL(a == b, true)
  TEST_EQUAL(a < b, true)
  TEST_EQUAL(b < a, false)

  // strict ordering: exactly one direction holds for these distinct masses
  TEST_EQUAL((a < b) != (b < a), true)
}
END_SECTION

START_SECTION(([EXTRA] std::hash<NuXLFragmentAdductDefinition> is consistent with operator==))
{
  EmpiricalFormula f("C5H10");
  std::hash<NuXLFragmentAdductDefinition> H;

  NuXLFragmentAdductDefinition a(f, "adductA", 70.0);
  NuXLFragmentAdductDefinition b(f, "adductA", 999.0); // == a (mass ignored)
  // equal objects must hash equally
  TEST_EQUAL(H(a) == H(b), true)

  // usable as an unordered_set key: a and b collapse to one entry, a different
  // formula adds a second
  std::unordered_set<NuXLFragmentAdductDefinition> s;
  s.insert(a);
  s.insert(b);
  TEST_EQUAL(s.size(), 1)
  s.insert(NuXLFragmentAdductDefinition(EmpiricalFormula("C6H12"), "adductA", 70.0));
  TEST_EQUAL(s.size(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
