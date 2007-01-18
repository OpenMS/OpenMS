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

#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/DATASTRUCTURES/HashMap.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ElementDB, "$Id$")

/////////////////////////////////////////////////////////////

const ElementDB* e_ptr = 0;
CHECK(static const ElementDB* getInstance())
	e_ptr = ElementDB::getInstance();
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK((const HashMap<String, const Element*>& getNames() const))
	HashMap<String, const Element*> names = e_ptr->getNames();
	const Element * e = e_ptr->getElement("Carbon");
	TEST_EQUAL(e, names["Carbon"])
	TEST_NOT_EQUAL(e, 0)
RESULT


CHECK((const HashMap<String, const Element*>& getSymbols() const))
	HashMap<String, const Element*> symbols = e_ptr->getSymbols();
	const Element * e = e_ptr->getElement("Carbon");
	TEST_EQUAL(e, symbols["C"])
	TEST_NOT_EQUAL(e, 0)
RESULT

CHECK((const HashMap<Size, const Element*>& getAtomicNumbers() const))
	HashMap<Size, const Element*> atomic_numbers = e_ptr->getAtomicNumbers();
	const Element * e = e_ptr->getElement("Carbon");
	TEST_EQUAL(e, atomic_numbers[6])
	TEST_NOT_EQUAL(e, 0)
RESULT

CHECK(const Element* getElement(const String& name) const)
	const Element * e1 = e_ptr->getElement("Hydrogen");
	const Element * e2 = e_ptr->getElement("H");
	TEST_EQUAL(e1, e2);
	TEST_NOT_EQUAL(e1, 0);
RESULT

CHECK(const Element* getElement(Size atomic_number) const)
	const Element * e1 = e_ptr->getElement("Carbon");
	const Element * e2 = e_ptr->getElement(6);
	TEST_EQUAL(e1, e2)
	TEST_NOT_EQUAL(e1, 0)
RESULT

CHECK(bool hasElement(const String& name) const)
	TEST_EQUAL(e_ptr->hasElement("Carbon"), true)
RESULT

CHECK(bool hasElement(Size atomic_number) const)
	TEST_EQUAL(e_ptr->hasElement(6), true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
