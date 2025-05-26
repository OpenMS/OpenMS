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
#include <OpenMS/CHEMISTRY/Isotope.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>


using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Element, "$Id$")

/////////////////////////////////////////////////////////////

Element* e_ptr = nullptr;
Element* e_nullPointer = nullptr;
START_SECTION(Element())
	e_ptr = new Element;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION(~Element())
	delete e_ptr;
END_SECTION

IsotopeDistribution dist;
string name("Name"), symbol("Symbol");
unsigned int atomic_number(43);
double average_weight(0.12345);
double mono_weight(0.123456789);

e_ptr = nullptr;
START_SECTION((Element(const string& name, const string& symbol, unsigned int atomic_number, double average_weight, double mono_weight, const IsotopeDistribution& isotopes)))
	e_ptr = new Element(name, symbol, atomic_number, average_weight, mono_weight, dist);	
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION(Element(const Element& element))
	Element copy(*e_ptr);
	TEST_EQUAL(*e_ptr == copy, true)
END_SECTION

delete e_ptr;
e_ptr = new Element;

START_SECTION(void setAtomicNumber(unsigned int atomic_number))
	e_ptr->setAtomicNumber(atomic_number);
	NOT_TESTABLE
END_SECTION

START_SECTION(UInt getAtomicNumber() const)
	TEST_EQUAL(e_ptr->getAtomicNumber(), atomic_number)
END_SECTION

START_SECTION(void setName(const string& name))
	e_ptr->setName(name);
	NOT_TESTABLE
END_SECTION

START_SECTION(const string& getName() const)
	TEST_EQUAL(e_ptr->getName(), name)
END_SECTION

START_SECTION(void setSymbol(const string& symbol))
	e_ptr->setSymbol(symbol);
	NOT_TESTABLE
END_SECTION

START_SECTION(const string& getSymbol() const)
	TEST_EQUAL(e_ptr->getSymbol(), symbol)
END_SECTION

START_SECTION(void setIsotopeDistribution(const IsotopeDistribution& isotopes))
	e_ptr->setIsotopeDistribution(dist);
	NOT_TESTABLE
END_SECTION

START_SECTION((const IsotopeDistribution& getIsotopeDistribution() const))
	TEST_EQUAL(e_ptr->getIsotopeDistribution() == dist, true)
END_SECTION

START_SECTION(void setAverageWeight(double weight))
	e_ptr->setAverageWeight(average_weight);
	NOT_TESTABLE
END_SECTION

START_SECTION(( void setIsotopes(const std::vector<const Isotope*>& isotopes) ))
  // Create isotopes with meaningful data to test cache functionality
  Isotope* iso1 = new Isotope("Carbon-12", "C", 6, 6, 12.0, 0.9893, -1, Isotope::DecayMode::NONE);
  Isotope* iso2 = new Isotope("Carbon-13", "C", 6, 7, 13.003355, 0.0107, -1, Isotope::DecayMode::NONE);
  std::vector<const Isotope*> isotopes = {iso1, iso2};
  
  // Store initial isotope distribution for comparison
  IsotopeDistribution initial_dist = e_ptr->getIsotopeDistribution();
  
  // Set isotopes and verify cache is updated
  e_ptr->setIsotopes(isotopes);
  
  // Verify the isotope distribution cache was updated correctly
  const IsotopeDistribution& updated_dist = e_ptr->getIsotopeDistribution();
  TEST_NOT_EQUAL(updated_dist == initial_dist, true) // Distribution should have changed
  TEST_EQUAL(updated_dist.size(), 2) // Should contain 2 isotopes
  
  // Verify the isotope distribution contains the correct masses and abundances
  auto container = updated_dist.getContainer();
  TEST_REAL_SIMILAR(container[0].getMZ(), 12.0)
  TEST_REAL_SIMILAR(container[0].getIntensity(), 0.9893)
  TEST_REAL_SIMILAR(container[1].getMZ(), 13.003355)
  TEST_REAL_SIMILAR(container[1].getIntensity(), 0.0107)
  
  // Clean up after the test
  for (auto iso : isotopes)
  {
    delete iso;
  }
END_SECTION

START_SECTION(( const std::vector<const Isotope*>& getIsotopes() const))
  TEST_EQUAL(e_ptr->getIsotopes().size(), 2)
END_SECTION

START_SECTION(double getAverageWeight() const)
	TEST_REAL_SIMILAR(e_ptr->getAverageWeight(), average_weight)
END_SECTION

START_SECTION(void setMonoWeight(double weight))
	e_ptr->setMonoWeight(2.333);
	NOT_TESTABLE
END_SECTION

START_SECTION(double getMonoWeight() const)
	TEST_REAL_SIMILAR(e_ptr->getMonoWeight(), 2.333)
END_SECTION

START_SECTION(virtual bool isIsotope())
  // Test with Element - should return false
  TEST_EQUAL( e_ptr->isIsotope(), false )
  
  // Test with actual Isotope - should return true
  Isotope isotope("Hydrogen-1", "H", 1, 0, 1.007825, 0.999885, -1, Isotope::DecayMode::NONE);
  TEST_EQUAL( isotope.isIsotope(), true )
END_SECTION

START_SECTION(Element& operator = (const Element& element))
	Element e = *e_ptr;
	TEST_EQUAL(e == *e_ptr, true)
END_SECTION

START_SECTION(bool operator != (const Element& element) const)
	Element e(*e_ptr);
	TEST_EQUAL(e != *e_ptr, false)
	e.setAverageWeight(0.54321);
	TEST_EQUAL(e != *e_ptr, true)
END_SECTION

START_SECTION(bool operator == (const Element& element) const)
	Element e(*e_ptr);
	TEST_EQUAL(e == *e_ptr, true)
	e.setAverageWeight(0.54321);
	TEST_EQUAL(e == *e_ptr, false)
END_SECTION

START_SECTION(bool operator < (const Element& element) const)
	const Element * h = ElementDB::getInstance()->getElement("H");
	const Element * c = ElementDB::getInstance()->getElement("Carbon");
	const Element * o = ElementDB::getInstance()->getElement("O");
	const Element * s = ElementDB::getInstance()->getElement("S");
	TEST_EQUAL(*h < *c, true)
	TEST_EQUAL(*c < *o, true)
	TEST_EQUAL(*c < *c, false)
	TEST_EQUAL(*s < *c, false)
END_SECTION


delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
