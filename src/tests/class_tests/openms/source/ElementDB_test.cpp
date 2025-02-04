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

ElementDB* e_ptr = nullptr;
ElementDB* e_nullPointer = nullptr;
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

START_SECTION(void addElement(const std::string& name,
                    const std::string& symbol,
                    const unsigned int an,
                    const std::map<unsigned int, double>& abundance,
                    const std::map<unsigned int, double>& mass,
                    bool replace_existing))
{
  const Element * oxygen = e_ptr->getElement(8);
  TEST_REAL_SIMILAR(oxygen->getAverageWeight(), 15.99940532316)
  map<unsigned int, double> oxygen_abundance = {{16u, 0.7}, {19u, 0.3}};
  map<unsigned int, double> oxygen_mass = {{16u, 15.994915000000001}, {19u, 19.01}};
  e_ptr->addElement("Oxygen", "O", 8u, oxygen_abundance, oxygen_mass, true); // true: replace existing

  const Element * new_oxygen = e_ptr->getElement(8);
  // ptr addresses cannot change, otherwise we are in trouble since EmpiricalFormula uses those
  TEST_EQUAL(oxygen, new_oxygen)
  TEST_REAL_SIMILAR(oxygen->getAverageWeight(), 16.8994405) // average weight has changed

  // cannot add invalid element (name and symbol conflict when compared existing element 
  // -- this would invalidate the lookup, since e_ptr->getSymbols().at("O")->getSymbol() == 'P'
  TEST_EXCEPTION(Exception::InvalidValue, e_ptr->addElement("Oxygen", "P", 8u, oxygen_abundance, oxygen_mass, true)) // true: replace existing

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
