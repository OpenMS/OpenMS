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
#include <OpenMS/CHEMISTRY/AdductInfo.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Constants.h>

using namespace OpenMS;
using namespace std;

START_TEST(AdductInfo, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

const double PROTON = Constants::PROTON_MASS_U;
TOLERANCE_ABSOLUTE(1e-3)

START_SECTION((static AdductInfo parseAdductString(const std::string& adduct)))
{
  AdductInfo a = AdductInfo::parseAdductString("M+H;1+");
  TEST_EQUAL(a.getCharge(), 1)
  TEST_EQUAL(a.getMolMultiplier(), 1)
  TEST_STRING_EQUAL(a.getName(), "M+H;1+")

  // multimer prefix
  AdductInfo dimer = AdductInfo::parseAdductString("2M+H;1+");
  TEST_EQUAL(dimer.getMolMultiplier(), 2)
  TEST_EQUAL(dimer.getCharge(), 1)

  // negative ion
  AdductInfo neg = AdductInfo::parseAdductString("M-H;1-");
  TEST_EQUAL(neg.getCharge(), -1)
  TEST_EQUAL(neg.getMolMultiplier(), 1)

  // metal adduct
  AdductInfo na = AdductInfo::parseAdductString("M+Na;1+");
  TEST_EQUAL(na.getCharge(), 1)
}
END_SECTION

START_SECTION((AdductInfo(const std::string& name, const EmpiricalFormula& adduct, int charge, UInt mol_multiplier = 1)))
{
  AdductInfo a("my_proton", EmpiricalFormula("H"), 1, 1);
  TEST_EQUAL(a.getCharge(), 1)
  TEST_STRING_EQUAL(a.getName(), "my_proton")
  TEST_EQUAL(a.getMolMultiplier(), 1)
  // built like [M+H]+ -> mass shift ~ 0 and m/z = M + proton
  TEST_REAL_SIMILAR(a.getMassShift(), 0.0)
  TEST_REAL_SIMILAR(a.getMZ(100.0), 100.0 + PROTON)

  // default mol_multiplier is 1
  AdductInfo b("na", EmpiricalFormula("Na"), 1);
  TEST_EQUAL(b.getMolMultiplier(), 1)
}
END_SECTION

START_SECTION((double getMZ(double neutral_mass) const))
{
  // [M+H]+   : m/z = M + proton
  TEST_REAL_SIMILAR(AdductInfo::parseAdductString("M+H;1+").getMZ(100.0), 100.0 + PROTON)
  // [2M+H]+  : m/z = 2M + proton
  TEST_REAL_SIMILAR(AdductInfo::parseAdductString("2M+H;1+").getMZ(100.0), 200.0 + PROTON)
  // [M-H]-   : m/z = M - proton
  TEST_REAL_SIMILAR(AdductInfo::parseAdductString("M-H;1-").getMZ(100.0), 100.0 - PROTON)
}
END_SECTION

START_SECTION((double getNeutralMass(double observed_mz) const))
{
  // getNeutralMass is the exact inverse of getMZ -> round-trips the input mass
  const double M = 523.25;
  for (const std::string& s : {std::string("M+H;1+"), std::string("M+Na;1+"),
                               std::string("2M+H;1+"), std::string("M-H;1-")})
  {
    AdductInfo a = AdductInfo::parseAdductString(s);
    TEST_REAL_SIMILAR(a.getNeutralMass(a.getMZ(M)), M)
  }
}
END_SECTION

START_SECTION((double getMassShift(bool use_avg_mass = false) const))
{
  // [M+H]+ : the added proton and the missing electron net out against the
  // added hydrogen -> mass shift ~ 0
  TEST_REAL_SIMILAR(AdductInfo::parseAdductString("M+H;1+").getMassShift(), 0.0)

  // [M+Na]+ : shift = mass(Na) - mass(H)
  const double expected = EmpiricalFormula("Na").getMonoWeight() - EmpiricalFormula("H").getMonoWeight();
  TEST_REAL_SIMILAR(AdductInfo::parseAdductString("M+Na;1+").getMassShift(), expected)

  // average-mass variant differs from the monoisotopic one for Na
  const double avg = AdductInfo::parseAdductString("M+Na;1+").getMassShift(true);
  TEST_EQUAL(avg > 0.0, true)
}
END_SECTION

START_SECTION((bool isCompatible(const EmpiricalFormula& db_entry) const))
{
  // [M-H]- removes an H -> only compatible with a candidate that contains H
  AdductInfo mh = AdductInfo::parseAdductString("M-H;1-");
  TEST_EQUAL(mh.isCompatible(EmpiricalFormula("C6H12O6")), true)
  TEST_EQUAL(mh.isCompatible(EmpiricalFormula("O2")), false)

  // a plus-only adduct subtracts nothing -> compatible with any compound
  AdductInfo na = AdductInfo::parseAdductString("M+Na;1+");
  TEST_EQUAL(na.isCompatible(EmpiricalFormula("O2")), true)
}
END_SECTION

START_SECTION((int getCharge() const))
{
  TEST_EQUAL(AdductInfo::parseAdductString("M+H;1+").getCharge(), 1)
  TEST_EQUAL(AdductInfo::parseAdductString("M-H;1-").getCharge(), -1)
}
END_SECTION

START_SECTION((const std::string& getName() const))
{
  TEST_STRING_EQUAL(AdductInfo::parseAdductString("M+Na;1+").getName(), "M+Na;1+")
}
END_SECTION

START_SECTION((const EmpiricalFormula& getEmpiricalFormula() const))
{
  // the stored adduct formula is the neutral adduct part (charge held separately)
  AdductInfo na = AdductInfo::parseAdductString("M+Na;1+");
  TEST_EQUAL(na.getEmpiricalFormula().getCharge(), 0)
}
END_SECTION

START_SECTION((UInt getMolMultiplier() const))
{
  TEST_EQUAL(AdductInfo::parseAdductString("M+H;1+").getMolMultiplier(), 1)
  TEST_EQUAL(AdductInfo::parseAdductString("2M+H;1+").getMolMultiplier(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
