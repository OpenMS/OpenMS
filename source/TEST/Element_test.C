// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/DATASTRUCTURES/String.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Element, "$Id$")

/////////////////////////////////////////////////////////////

Element* e_ptr = 0;
CHECK(Element())
	e_ptr = new Element;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~Element())
	delete e_ptr;
RESULT

IsotopeDistribution dist;
String name("Name"), symbol("Symbol");
Size atomic_number(43);
Real average_weight(0.12345);
Real mono_weight(0.123456789);

e_ptr = 0;
CHECK((Element(const String& name, const String& symbol, Size atomic_number, Real average_weight, Real mono_weight, const IsotopeDistribution& isotopes)))
	e_ptr = new Element(name, symbol, atomic_number, average_weight, mono_weight, dist);	
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(Element(const Element& element))
	Element copy(*e_ptr);
	TEST_EQUAL(*e_ptr == copy, true)
RESULT

delete e_ptr;
e_ptr = new Element;

CHECK(void setAtomicNumber(Size atomic_number))
	e_ptr->setAtomicNumber(atomic_number);
RESULT

CHECK(Size getAtomicNumber() const)
	TEST_EQUAL(e_ptr->getAtomicNumber(), atomic_number)
RESULT

CHECK(void setName(const String& name))
	e_ptr->setName(name);
RESULT

CHECK(const String& getName() const)
	TEST_EQUAL(e_ptr->getName(), name)
RESULT

CHECK(void setSymbol(const String& symbol))
	e_ptr->setSymbol(symbol);
RESULT

CHECK(const String& getSymbol() const)
	TEST_EQUAL(e_ptr->getSymbol(), symbol)
RESULT

CHECK(void setIsotopeDistribution(const IsotopeDistribution& isotopes))
	e_ptr->setIsotopeDistribution(dist);
RESULT

CHECK((const IsotopeDistribution& getIsotopeDistribution() const))
	TEST_EQUAL(e_ptr->getIsotopeDistribution() == dist, true)
RESULT

CHECK(void setAverageWeight(Real weight))
	e_ptr->setAverageWeight(average_weight);
RESULT

CHECK(Real getAverageWeight() const)
	TEST_REAL_EQUAL(e_ptr->getAverageWeight(), average_weight)
RESULT

CHECK(void setMonoWeight(Real weight))
	e_ptr->setMonoWeight(2.333);
RESULT

CHECK(Real getMonoWeight() const)
	TEST_REAL_EQUAL(e_ptr->getMonoWeight(), 2.333)
RESULT

CHECK(Element& operator = (const Element& element))
	Element e = *e_ptr;
	TEST_EQUAL(e == *e_ptr, true)
RESULT

CHECK(bool operator != (const Element& element) const)
	Element e(*e_ptr);
	TEST_EQUAL(e != *e_ptr, false)
	e.setAverageWeight(0.54321);
	TEST_EQUAL(e != *e_ptr, true)
RESULT

CHECK(bool operator == (const Element& element) const)
	Element e(*e_ptr);
	TEST_EQUAL(e == *e_ptr, true)
	e.setAverageWeight(0.54321);
	TEST_EQUAL(e == *e_ptr, false)
RESULT

CHECK(friend std::ostream& operator << (std::ostream& os, const Element& element))
	// ???
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
