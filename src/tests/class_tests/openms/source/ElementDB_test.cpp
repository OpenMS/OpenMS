// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
  const Isotope * e1 = e_ptr->getIsotope("(238)U");
  const Isotope * e2 = e_ptr->getIsotope("(238)Uranium");
  TEST_EQUAL(e1, e2);
  TEST_NOT_EQUAL(e1, elem_nullPointer);
  TEST_EQUAL(e1->getNeutrons(), 146);
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

START_SECTION([extra] output generation)
{
  std::cout << *(e_ptr->getElement(8)) << std::endl;
  std::cout << *(e_ptr->getElement(6)) << std::endl;
  std::cout << *(e_ptr->getElement(92)) << std::endl;
  NOT_TESTABLE
}
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

  TEST_EQUAL(e_ptr->getElement(800), elem_nullPointer)
  e_ptr->addElement("NewElement", "NE", 800u, oxygen_abundance, oxygen_mass, false);
  const Element * new_ele = e_ptr->getElement(800);
  TEST_REAL_SIMILAR(new_ele->getAverageWeight(), 16.8994405) // average weight of new element

  // Test that we cannot add elements twice with replace=false 
  TEST_EXCEPTION(Exception::IllegalArgument, e_ptr->addElement("NewElement", "NE", 800u, oxygen_abundance, oxygen_mass, false))

  // Test that we can add elements twice with replace=true 
  e_ptr->addElement("NewElement", "NE", 800u, oxygen_abundance, oxygen_mass, true);

  // cannot add invalid element (name and symbol conflict when compared existing element 
  // -- this would invalidate the lookup, since e_ptr->getSymbols().at("O")->getSymbol() == 'P'
  TEST_EXCEPTION(Exception::InvalidValue, e_ptr->addElement("Oxygen", "P", 8u, oxygen_abundance, oxygen_mass, true)) // true: replace existing
}
END_SECTION

START_SECTION( void addIsotope(const std::string& name, const std::string& symbol, const unsigned int an, double abundance, double mass, double half_life, Isotope::DecayMode decay, bool replace_existing))
{
  const Isotope * iso1 = e_ptr->getIsotope("(238)U");
  TEST_REAL_SIMILAR(iso1->getAbundance(), 0.992742) // test natural abundance 
  e_ptr->addIsotope("Uranium", "U", 92u, 1.3, 238.05, 1e5, Isotope::DecayMode::UNKNOWN, true);

  const Isotope * iso2 = e_ptr->getIsotope("(238)U");
  // ptr addresses cannot change, otherwise we are in trouble since EmpiricalFormula uses those
  TEST_EQUAL(iso1, iso2)
  TEST_REAL_SIMILAR(iso1->getAbundance(), 1.3) // natural abundance has changed

  // we have now managed to have 130% natural abundance for Uranium
  // NOTE: this is a major problem for average weight calculations etc
  const Element* element = e_ptr->getElement(92);
  double sum = 0;
  for (auto& iso : element->getIsotopeDistribution()) {sum += iso.getIntensity();}
  TEST_REAL_SIMILAR(sum, 1.30725795222315);

  // new Uranium isotope added
  TEST_EQUAL(e_ptr->getIsotope("(314)C") == nullptr, true)
  int nr_isotopes = element->getIsotopes().size();
  e_ptr->addIsotope("Uranium", "U", 92u, 0.3, 314.0, 1e5, Isotope::DecayMode::UNKNOWN, false);
  const Isotope * new_iso = e_ptr->getIsotope("(314)U");
  TEST_EQUAL( new_iso != nullptr, true)
  TEST_EQUAL( new_iso->getElement(), element)
  TEST_EQUAL( element->getIsotopes().size(), nr_isotopes+1) // increased number of isotopes by one
  TEST_EQUAL( element->getIsotopeDistribution().getContainer().size(), nr_isotopes+1) // increased number of isotopes by one

  // we have now managed to have 160% natural abundance for Uranium
  // NOTE: this is a major problem for average weight calculations etc
  sum = 0;
  for (auto& iso : element->getIsotopeDistribution()) {sum += iso.getIntensity();}
  TEST_REAL_SIMILAR(sum, 1.60725795222315);

  // Test that we cannot add isotopes for elements that dont exist
  TEST_EXCEPTION(Exception::IllegalArgument, e_ptr->addIsotope("NewElement", "NE", 300, 0.5, 404, 100, Isotope::DecayMode::UNKNOWN, false))

  {
    // Test that we cannot add twice with replace=false
    e_ptr->addIsotope("NewElement", "NE", 800, 0.5, 404, 100, Isotope::DecayMode::UNKNOWN, false);
    TEST_EXCEPTION(Exception::IllegalArgument, e_ptr->addIsotope("NewElement", "NE", 800, 0.5, 404, 100, Isotope::DecayMode::UNKNOWN, false))
    const Isotope* new_iso = e_ptr->getIsotope("(404)NE");
    TEST_EQUAL( new_iso != nullptr, true);
    TEST_REAL_SIMILAR( new_iso->getAbundance(), 0.5)

    // we should be able to add the same one again with replace=true
    e_ptr->addIsotope("NewElement", "NE", 800, 0.6, 404, 100, Isotope::DecayMode::UNKNOWN, true);
    new_iso = e_ptr->getIsotope("(404)NE");
    TEST_EQUAL( new_iso != nullptr, true);
    TEST_REAL_SIMILAR( new_iso->getAbundance(), 0.6)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
