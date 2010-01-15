// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Andreas Bertsch $
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

START_TEST(Element, "$Id: Element_test.C 5908 2009-08-26 13:44:26Z marc_sturm $")

/////////////////////////////////////////////////////////////

Element* e_ptr = 0;
START_SECTION(Element())
	e_ptr = new Element;
	TEST_NOT_EQUAL(e_ptr, 0)
END_SECTION

START_SECTION(~Element())
	delete e_ptr;
END_SECTION

IsotopeDistribution dist;
String name("Name"), symbol("Symbol");
UInt atomic_number(43);
DoubleReal average_weight(0.12345);
DoubleReal mono_weight(0.123456789);

e_ptr = 0;
START_SECTION((Element(const String& name, const String& symbol, UInt atomic_number, DoubleReal average_weight, DoubleReal mono_weight, const IsotopeDistribution& isotopes)))
	e_ptr = new Element(name, symbol, atomic_number, average_weight, mono_weight, dist);	
	TEST_NOT_EQUAL(e_ptr, 0)
END_SECTION

START_SECTION(Element(const Element& element))
	Element copy(*e_ptr);
	TEST_EQUAL(*e_ptr == copy, true)
END_SECTION

delete e_ptr;
e_ptr = new Element;

START_SECTION(void setAtomicNumber(UInt atomic_number))
	e_ptr->setAtomicNumber(atomic_number);
	NOT_TESTABLE
END_SECTION

START_SECTION(UInt getAtomicNumber() const)
	TEST_EQUAL(e_ptr->getAtomicNumber(), atomic_number)
END_SECTION

START_SECTION(void setName(const String& name))
	e_ptr->setName(name);
	NOT_TESTABLE
END_SECTION

START_SECTION(const String& getName() const)
	TEST_EQUAL(e_ptr->getName(), name)
END_SECTION

START_SECTION(void setSymbol(const String& symbol))
	e_ptr->setSymbol(symbol);
	NOT_TESTABLE
END_SECTION

START_SECTION(const String& getSymbol() const)
	TEST_EQUAL(e_ptr->getSymbol(), symbol)
END_SECTION

START_SECTION(void setIsotopeDistribution(const IsotopeDistribution& isotopes))
	e_ptr->setIsotopeDistribution(dist);
	NOT_TESTABLE
END_SECTION

START_SECTION((const IsotopeDistribution& getIsotopeDistribution() const))
	TEST_EQUAL(e_ptr->getIsotopeDistribution() == dist, true)
END_SECTION

START_SECTION(void setAverageWeight(DoubleReal weight))
	e_ptr->setAverageWeight(average_weight);
	NOT_TESTABLE
END_SECTION

START_SECTION(DoubleReal getAverageWeight() const)
	TEST_REAL_SIMILAR(e_ptr->getAverageWeight(), average_weight)
END_SECTION

START_SECTION(void setMonoWeight(DoubleReal weight))
	e_ptr->setMonoWeight(2.333);
	NOT_TESTABLE
END_SECTION

START_SECTION(DoubleReal getMonoWeight() const)
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

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
