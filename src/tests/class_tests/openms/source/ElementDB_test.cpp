// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ElementDB, "$Id$")

/////////////////////////////////////////////////////////////

const ElementDB* e_ptr = nullptr;
const ElementDB* e_nullPointer = nullptr;
const Element * elem_nullPointer = nullptr;

START_SECTION([EXTRA] multithreaded example)
{
  int nr_iterations (100);
  int test = 0;
#pragma omp parallel for reduction(+: test)
  for (int k = 1; k < nr_iterations + 1; k++)
  {
    auto edb = ElementDB::getInstance();
    const Element * e1 = edb->getElement("Carbon");
    test += e1->getAtomicNumber();
  }
  TEST_EQUAL(test, 6 * 100)
}
END_SECTION

START_SECTION(static const ElementDB* getInstance())
  e_ptr = ElementDB::getInstance();
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((const unordered_map<string, const Element*>& getNames() const))
  unordered_map<string, const Element*> names = e_ptr->getNames();
  const Element * e = e_ptr->getElement("Carbon");
  TEST_EQUAL(e, names["Carbon"])
  TEST_NOT_EQUAL(e, elem_nullPointer)
END_SECTION


START_SECTION((const unordered_map<string, const Element*>& getSymbols() const))
  unordered_map<string, const Element*> symbols = e_ptr->getSymbols();
  const Element * e = e_ptr->getElement("Carbon");
  TEST_EQUAL(e, symbols["C"])
  TEST_NOT_EQUAL(e, elem_nullPointer)
END_SECTION

START_SECTION((const unordered_map<unsigned int, const Element*>& getAtomicNumbers() const))
  unordered_map<unsigned int, const Element*> atomic_numbers = e_ptr->getAtomicNumbers();
  const Element * e = e_ptr->getElement("Carbon");
  TEST_EQUAL(e, atomic_numbers[6])
  TEST_NOT_EQUAL(e, elem_nullPointer)
END_SECTION

START_SECTION(const Element* getElement(const string& name) const)
  const Element * e1 = e_ptr->getElement("Hydrogen");
  const Element * e2 = e_ptr->getElement("H");
  TEST_EQUAL(e1, e2);
  TEST_NOT_EQUAL(e1, elem_nullPointer);
END_SECTION

START_SECTION(const Element* getElement(unsigned int atomic_number) const)
  const Element * e1 = e_ptr->getElement("Carbon");
  const Element * e2 = e_ptr->getElement(6);
  TEST_EQUAL(e1, e2)
  TEST_NOT_EQUAL(e1, elem_nullPointer)
END_SECTION

START_SECTION(bool hasElement(const string& name) const)
  TEST_EQUAL(e_ptr->hasElement("Carbon"), true)
END_SECTION

START_SECTION(bool hasElement(unsigned int atomic_number) const)
  TEST_EQUAL(e_ptr->hasElement(6), true)
END_SECTION

START_SECTION([EXTRA] iridium has its own isotopes)
{
  // natural iridium: 191Ir 190.960591 u (37.23 %), 193Ir 192.962924 u (62.77 %)
  const Element* ir = e_ptr->getElement("Ir");
  TEST_EQUAL(ir->getAtomicNumber(), 77)
  const IsotopeDistribution& iso = ir->getIsotopeDistribution();
  TEST_EQUAL(iso.size(), 2)
  TEST_REAL_SIMILAR(iso[0].getMZ(), 190.960591)
  TEST_REAL_SIMILAR(iso[0].getIntensity(), 0.3723)
  TEST_REAL_SIMILAR(iso[1].getMZ(), 192.962924)
  TEST_REAL_SIMILAR(iso[1].getIntensity(), 0.6277)
  TEST_REAL_SIMILAR(ir->getMonoWeight(), 192.962924) // most abundant isotope
  TEST_REAL_SIMILAR(ir->getAverageWeight(), 192.2174554) // 0.3723 * 190.960591 + 0.6277 * 192.962924
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
