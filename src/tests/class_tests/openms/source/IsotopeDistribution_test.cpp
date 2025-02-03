// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//

///////////////////////////
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
///////////////////////////

// More headers

#include <iostream>
#include <iterator>
#include <utility>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

/////////////////////////////////////////////////////////////

START_TEST(IsotopeDistribution, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

using namespace OpenMS;
using namespace std;

IsotopeDistribution* nullPointer = nullptr;

START_SECTION(CoarseIsotopePatternGenerator())
  IsotopeDistribution* ptr = nullptr;
  ptr = new IsotopeDistribution();
  Size container_size = ptr->size();
  TEST_EQUAL(container_size, 1)
  TEST_NOT_EQUAL(ptr, nullPointer)
  delete ptr;

  // Ensure that IsotopeDistribution has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(IsotopeDistribution(std::declval<IsotopeDistribution&&>())), true)
END_SECTION

IsotopeDistribution* iso = new IsotopeDistribution();

START_SECTION(IsotopeDistribution(const IsotopeDistribution& isotope_distribution))
  IsotopeDistribution copy;
  copy = *iso;
  for (Size i = 0; i != copy.getContainer().size(); ++i)
  {
    TEST_EQUAL(copy.getContainer()[i].getMZ(), iso->getContainer()[i].getMZ())
    TEST_EQUAL(copy.getContainer()[i].getIntensity(), iso->getContainer()[i].getIntensity())
  }
  TEST_EQUAL(copy.getMin(), iso->getMin())
  TEST_EQUAL(copy.getMax(), iso->getMax())
  TEST_EQUAL(copy.size(), iso->size())
END_SECTION

START_SECTION(~IsotopeDistribution())
  IsotopeDistribution* ptr = new IsotopeDistribution();
  delete ptr;
END_SECTION

START_SECTION(IsotopeDistribution& operator = (const CoarseIsotopePatternGenerator& isotope_distribution))
  IsotopeDistribution copy;
  copy = *iso;
  for (Size i = 0; i != copy.getContainer().size(); ++i)
  {
    TEST_EQUAL(copy.getContainer()[i].getMZ(), iso->getContainer()[i].getMZ())
    TEST_EQUAL(copy.getContainer()[i].getIntensity(), iso->getContainer()[i].getIntensity())
  }
  TEST_EQUAL(copy.getMin(), iso->getMin())
  TEST_EQUAL(copy.getMax(), iso->getMax())
  TEST_EQUAL(copy.size(), iso->size())
END_SECTION

START_SECTION(IsotopeDistribution& operator < (const CoarseIsotopePatternGenerator& isotope_distribution))
  IsotopeDistribution iso1, iso2;
  TEST_EQUAL(iso1 < iso2, false)
  IsotopeDistribution iso3(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11))),
    iso4(EmpiricalFormula("C5").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
  TEST_EQUAL(iso3 < iso4, true)

  IsotopeDistribution iso5(EmpiricalFormula("C5").getIsotopeDistribution(CoarseIsotopePatternGenerator(1)));
  IsotopeDistribution iso6(EmpiricalFormula("C5").getIsotopeDistribution(CoarseIsotopePatternGenerator(1000)));
  TEST_EQUAL(iso5 < iso6, true)

  IsotopeDistribution iso7(EmpiricalFormula("C5").getIsotopeDistribution(CoarseIsotopePatternGenerator(11, true)));
  IsotopeDistribution iso8(EmpiricalFormula("C5").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
  // iso7 should be less because its second isotope's mass is 61 (atomic number), while for iso8 it is 61.003 (expected mass)
  TEST_EQUAL(iso7 < iso8, true)
END_SECTION


START_SECTION(bool operator==(const IsotopeDistribution &isotope_distribution) const)
  IsotopeDistribution iso1, iso2;
  TEST_TRUE(iso1 == iso2)
  IsotopeDistribution iso3(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11))),
    iso4(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
  TEST_TRUE(iso3 == iso4)

  IsotopeDistribution iso5(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11, true))),
    iso6(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
  // the masses should be different
  TEST_EQUAL(iso5 == iso6, false)
END_SECTION

START_SECTION(void set(const ContainerType &distribution))
  IsotopeDistribution iso1(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11))), iso2;
  TEST_EQUAL(iso1 == iso2, false)
  IsotopeDistribution::ContainerType container = iso1.getContainer();
  iso2.set(container);
  TEST_EQUAL(iso1.getContainer() == iso2.getContainer(), true)
  TEST_TRUE(iso1 == iso2)
END_SECTION

START_SECTION(const ContainerType& getContainer() const)
  NOT_TESTABLE
END_SECTION

START_SECTION(Size getMax() const)
  IsotopeDistribution iso(EmpiricalFormula("H2").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
  TEST_REAL_SIMILAR(iso.getMax(), 6.02907)
  IsotopeDistribution iso2(EmpiricalFormula("H2").getIsotopeDistribution(CoarseIsotopePatternGenerator(11, true)));
  TEST_EQUAL(iso2.getMax(), 6)

  iso.insert(11.2, 2.0);
  iso.insert(10.2, 2.0);
  TEST_REAL_SIMILAR(iso.getMax(), 11.2)
END_SECTION

START_SECTION(Size getMin() const)
  IsotopeDistribution iso(EmpiricalFormula("H2").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
  TEST_REAL_SIMILAR(iso.getMin(), 2.01565)
  IsotopeDistribution iso2(EmpiricalFormula("H2").getIsotopeDistribution(CoarseIsotopePatternGenerator(11, true)));
  TEST_EQUAL(iso2.getMin(), 2)
  IsotopeDistribution iso3(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
  TEST_REAL_SIMILAR(iso3.getMin(), 48)
  IsotopeDistribution iso4(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11, true)));
  TEST_EQUAL(iso4.getMin(), 48)

  iso.insert(1.2, 2.0);
  iso.insert(10.2, 2.0);
  TEST_REAL_SIMILAR(iso.getMin(), 1.2)
END_SECTION

START_SECTION(Size getMostAbundant() const)
  IsotopeDistribution iso(EmpiricalFormula("C1").getIsotopeDistribution(CoarseIsotopePatternGenerator(11, true)));
  // The most abundant isotope is the monoisotope
  TEST_EQUAL(iso.getMostAbundant().getMZ(), 12)
  IsotopeDistribution iso2(EmpiricalFormula("C100").getIsotopeDistribution(CoarseIsotopePatternGenerator(11, true)));
  // In this case, the most abundant isotope isn't the monoisotope
  TEST_EQUAL(iso2.getMostAbundant().getMZ(), 1201)
  // Empty distribution
  iso2.clear();
  TEST_EQUAL(iso2.getMostAbundant().getMZ(), 0);
  TEST_EQUAL(iso2.getMostAbundant().getIntensity(), 1);
END_SECTION

START_SECTION(Size size() const)
  IsotopeDistribution iso1, iso2(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
  TEST_EQUAL(iso1.size(), 1)
  TEST_EQUAL(iso2.size(), 5)
END_SECTION

START_SECTION(void clear())
  IsotopeDistribution iso2(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
  TEST_EQUAL(iso2.size(), 5)
  iso2.clear();
  TEST_EQUAL(iso2.size(), 0)
END_SECTION


START_SECTION(void trimRight(double cutoff))
  IsotopeDistribution iso(EmpiricalFormula("C160").getIsotopeDistribution(CoarseIsotopePatternGenerator(10)));
  TEST_NOT_EQUAL(iso.size(),3)
  iso.trimRight(0.2);
  TEST_EQUAL(iso.size(),3)
END_SECTION

START_SECTION(void trimLeft(double cutoff))
  IsotopeDistribution iso(EmpiricalFormula("C160").getIsotopeDistribution(CoarseIsotopePatternGenerator(10)));
  iso.trimRight(0.2);
  iso.trimLeft(0.2);
  TEST_EQUAL(iso.size(),2)
END_SECTION

START_SECTION(void renormalize())
  IsotopeDistribution iso(EmpiricalFormula("C160").getIsotopeDistribution(CoarseIsotopePatternGenerator(10)));
  iso.trimRight(0.2);
  iso.trimLeft(0.2);
  iso.renormalize();
  double sum = 0;
  for (IsotopeDistribution::ConstIterator it = iso.begin(); it != iso.end(); ++it)
  {
    sum += it->getIntensity();
  }

  TEST_REAL_SIMILAR(sum, 1.0)
END_SECTION



START_SECTION(bool operator!=(const IsotopeDistribution &isotope_distribution) const)
  IsotopeDistribution iso1;
  IsotopeDistribution iso2;
  TEST_EQUAL(iso1 != iso2, false)
  IsotopeDistribution iso3(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11))),
                      iso4(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
  TEST_EQUAL(iso3 != iso4, false)
  TEST_FALSE(iso2 == iso3)
END_SECTION

START_SECTION(Iterator begin())
  NOT_TESTABLE
END_SECTION

START_SECTION(Iterator end())
  NOT_TESTABLE
END_SECTION

START_SECTION(ConstIterator begin() const)
  NOT_TESTABLE
END_SECTION

START_SECTION(ConstIterator end() const)
  NOT_TESTABLE
END_SECTION

START_SECTION(ReverseIterator rbegin())
  NOT_TESTABLE
END_SECTION

START_SECTION(ReverseIterator rend())
  NOT_TESTABLE
END_SECTION

START_SECTION(ConstReverseIterator rbegin() const)
  NOT_TESTABLE
END_SECTION

START_SECTION(ConstReverseIterator rend() const)
  NOT_TESTABLE
END_SECTION

delete iso;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#pragma clang diagnostic pop

