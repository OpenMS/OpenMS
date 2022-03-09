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

START_SECTION(const Isotope* getIsotope(const string& name) const)
  const Isotope * e1 = e_ptr->getIsotope("(14)C");
  const Isotope * e2 = e_ptr->getIsotope("(14)Carbon");
  TEST_EQUAL(e1, e2);
  TEST_NOT_EQUAL(e1, elem_nullPointer);
  TEST_EQUAL(e1->getNeutrons(), 8);
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


START_SECTION(void addElement(const std::string& name, const std::string& symbol, const unsigned int an, const std::map<unsigned int, double>& abundance, const std::map<unsigned int, double>& mass, bool replace_existing))
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

  TEST_EQUAL(e_ptr->getElement(800), nullptr)
  e_ptr->addElement("NewElement", "NE", 800u, oxygen_abundance, oxygen_mass, false);
  const Element * new_ele = e_ptr->getElement(800);
  TEST_REAL_SIMILAR(new_ele->getAverageWeight(), 16.8994405) // average weight of new element
}
END_SECTION

START_SECTION( void addIsotope(const std::string& name, const std::string& symbol, const unsigned int an, double abundance, double mass, double half_life, Isotope::DecayMode decay, bool replace_existing))
{
  const Isotope * carbon14 = e_ptr->getIsotope("(14)C");
  e_ptr->addIsotope("Carbon", "C", 6u, 0.3, 14.0, 1e5, Isotope::DecayMode::UNKNOWN, true);

  const Isotope * new_carbon14 = e_ptr->getIsotope("(14)C");
  // ptr addresses cannot change, otherwise we are in trouble since EmpiricalFormula uses those
  TEST_EQUAL(carbon14, new_carbon14)
  TEST_REAL_SIMILAR(carbon14->getAbundance(), 0.3) // natural abundance has changed

  // we have now managed to have 130% natural abundance for Carbon
  // NOTE: this is a major problem for average weight calculations etc
  const Element* carbon = e_ptr->getElement(6);
  double sum = 0;
  for (auto& iso : carbon->getIsotopeDistribution()) {sum += iso.getIntensity();}
  TEST_REAL_SIMILAR(sum, 1.30000002495944);

  // new carbon isotope added
  TEST_EQUAL(e_ptr->getIsotope("(114)C"), nullptr)
  int nr_isotopes = carbon->getIsotopes().size();
  e_ptr->addIsotope("Carbon", "C", 6u, 0.3, 114.0, 1e5, Isotope::DecayMode::UNKNOWN, false);
  const Isotope * new_iso = e_ptr->getIsotope("(114)C");
  TEST_EQUAL( new_iso->getElement(), carbon)
  TEST_EQUAL( carbon->getIsotopes().size(), nr_isotopes+1)

  // we have now managed to have 160% natural abundance for Carbon
  // NOTE: this is a major problem for average weight calculations etc
  sum = 0;
  for (auto& iso : carbon->getIsotopeDistribution()) {sum += iso.getIntensity();}
  TEST_REAL_SIMILAR(sum, 1.60000002495944);
}
END_SECTION

START_SECTION([extra] cout)
{
  std::cout << *(e_ptr->getElement(8)) << std::endl;
  std::cout << *(e_ptr->getElement(6)) << std::endl;
  std::cout << *(e_ptr->getElement(92)) << std::endl;
  NOT_TESTABLE
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
