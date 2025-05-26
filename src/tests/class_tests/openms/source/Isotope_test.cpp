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

#include <OpenMS/CHEMISTRY/Isotope.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Isotope, "$Id$")

/////////////////////////////////////////////////////////////

Isotope* e_ptr = nullptr;
Isotope* e_nullPointer = nullptr;
START_SECTION(Isotope())
	e_ptr = new Isotope;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION(~Isotope())
	delete e_ptr;
END_SECTION

IsotopeDistribution dist;
string name("Name"), symbol("Symbol");
unsigned int atomic_number(15);
double mono_weight(0.123456789);

e_ptr = nullptr;
START_SECTION((Isotope(const std::string & name,
            const std::string & symbol,
            unsigned int atomic_number,
            unsigned int neutrons,
            double mono_weight,
            double abundance,
            double half_life,
            Isotope::DecayMode dm)))
	e_ptr = new Isotope(name, symbol, atomic_number, 10u, mono_weight, 0.6, 42, Isotope::DecayMode::UNKNOWN);	
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION( const Element* getElement() const)
  // Create a specific imaginary isotope of Phosphorus for this test
  Isotope test_phosphorus("Phosphorus-32", "P", 15, 17, 31.973907, 0.0, 1233360, Isotope::DecayMode::BETA_MINUS);
  
  const Element* el = test_phosphorus.getElement();
  TEST_NOT_EQUAL(el, e_nullPointer)
  TEST_EQUAL(el->getSymbol(), "P")
  TEST_EQUAL(el->getAtomicNumber(), 15)
  TEST_EQUAL(el->getName(), "Phosphorus-32")
  
  // Also test with the existing e_ptr for completeness
  const Element* el_existing = e_ptr->getElement();
  TEST_NOT_EQUAL(el_existing, e_nullPointer)
  TEST_EQUAL(el_existing->getSymbol(), "Symbol")  // This matches what was actually set
END_SECTION

START_SECTION(Isotope(const Isotope& Isotope))
	Isotope copy(*e_ptr);
	TEST_EQUAL(*e_ptr == copy, true)
END_SECTION

delete e_ptr;
e_ptr = new Isotope;

START_SECTION(void setHalfLife(double hl))
	e_ptr->setHalfLife(8.6);
	NOT_TESTABLE
END_SECTION

START_SECTION(double getHalfLife() const)
	TEST_REAL_SIMILAR(e_ptr->getHalfLife(), 8.6)
END_SECTION

START_SECTION(void setAbundance(double hl))
	e_ptr->setAbundance(0.6);
	NOT_TESTABLE
END_SECTION

START_SECTION(double getAbundance() const)
	TEST_REAL_SIMILAR(e_ptr->getAbundance(), 0.6)
END_SECTION

START_SECTION(void setNeutrons(int hl))
	e_ptr->setNeutrons(10);
	NOT_TESTABLE
END_SECTION

START_SECTION(int getNeutrons() const)
	TEST_EQUAL(e_ptr->getNeutrons(), 10)
END_SECTION

START_SECTION(void setDecayMode(int hl))
	e_ptr->setDecayMode(Isotope::DecayMode::ALPHA);
	NOT_TESTABLE
END_SECTION

START_SECTION(int getDecayMode() const)
	TEST_EQUAL(e_ptr->getDecayMode(), Isotope::DecayMode::ALPHA)
END_SECTION

START_SECTION( virtual bool isIsotope() )
  TEST_EQUAL(e_ptr->isIsotope(), true)
END_SECTION

START_SECTION(bool isStable() const)
  TEST_EQUAL(e_ptr->isStable(), false)
END_SECTION

START_SECTION(Isotope& operator = (const Isotope& Isotope))
	Isotope e;
	e = *e_ptr;
	TEST_EQUAL(e == *e_ptr, true)
END_SECTION

START_SECTION(bool operator != (const Isotope& Isotope) const)
	Isotope e(*e_ptr);
	TEST_EQUAL(e != *e_ptr, false)
	e.setAverageWeight(0.54321);
	TEST_EQUAL(e != *e_ptr, true)
END_SECTION

START_SECTION(bool operator == (const Isotope& Isotope) const)
	Isotope e(*e_ptr);
	TEST_EQUAL(e == *e_ptr, true)
	e.setAverageWeight(0.54321);
	TEST_EQUAL(e == *e_ptr, false)
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


