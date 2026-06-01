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
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>

#include <unordered_set>
#include <unordered_map>


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
// Hash tests
/////////////////////////////////////////////////////////////

START_SECTION(([EXTRA] std::hash<Element>))
{
  // Test that equal objects have equal hashes
  IsotopeDistribution dist1;
  dist1.insert(12.0, 0.989);
  dist1.insert(13.003355, 0.011);

  Element e1("Carbon", "C", 6, 12.0107, 12.0, dist1);
  Element e2("Carbon", "C", 6, 12.0107, 12.0, dist1);

  std::hash<Element> hasher;
  TEST_EQUAL(e1 == e2, true)
  TEST_EQUAL(hasher(e1), hasher(e2))

  // Test that different objects (likely) have different hashes
  Element e3("Nitrogen", "N", 7, 14.0067, 14.003074, IsotopeDistribution());
  TEST_EQUAL(e1 == e3, false)
  // Note: Different objects should have different hashes (with very high probability)
  TEST_NOT_EQUAL(hasher(e1), hasher(e3))

  // Test consistency: same object hashes to same value
  TEST_EQUAL(hasher(e1), hasher(e1))

  // Test use in unordered_set
  std::unordered_set<Element> element_set;
  element_set.insert(e1);
  element_set.insert(e2); // Should not increase size (duplicate)
  element_set.insert(e3);
  TEST_EQUAL(element_set.size(), 2)
  TEST_EQUAL(element_set.count(e1), 1)
  TEST_EQUAL(element_set.count(e3), 1)

  // Test use in unordered_map
  std::unordered_map<Element, std::string> element_map;
  element_map[e1] = "first carbon";
  element_map[e2] = "second carbon"; // Should overwrite
  element_map[e3] = "nitrogen";
  TEST_EQUAL(element_map.size(), 2)
  TEST_EQUAL(element_map[e1], "second carbon")
  TEST_EQUAL(element_map[e3], "nitrogen")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
