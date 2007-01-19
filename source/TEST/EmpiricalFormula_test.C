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

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ElementDB, "$Id$")

/////////////////////////////////////////////////////////////

EmpiricalFormula* e_ptr = 0;
CHECK(EmpiricalFormula())
	e_ptr = new EmpiricalFormula;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~EmpiricalFormula())
	delete e_ptr;
RESULT

CHECK(EmpiricalFormula(const String& rhs) throw(Exception::ParseError))
	e_ptr = new EmpiricalFormula("C4");
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(EmpiricalFormula(const EmpiricalFormula& rhs))
	EmpiricalFormula ef(*e_ptr);
	TEST_EQUAL(ef == *e_ptr, true)
RESULT

CHECK(EmpiricalFormula(Size number, const Element* element, SignedInt charge = 0))
	EmpiricalFormula ef(4, e_ptr->getElement("C"));
	TEST_EQUAL(ef == *e_ptr, true)
	TEST_EQUAL(ef.getCharge(), 0)
RESULT

CHECK(const Element* getElement(Size atomic_number) const)
	const Element* e = e_ptr->getElement(6);
	TEST_EQUAL(e->getSymbol(), "C")
RESULT

CHECK(const Element* getElement(const String& name) const)
	const Element* e = e_ptr->getElement("C");
	TEST_EQUAL(e->getSymbol(), "C")
RESULT

CHECK(Size getNumberOf(Size atomic_number) const)
	Size num1 = e_ptr->getNumberOf(6);
	TEST_EQUAL(num1, 4);
RESULT

CHECK(Size getNumberOf(const String& name) const)
	Size num2 = e_ptr->getNumberOf("C");
	TEST_EQUAL(num2, 4);
RESULT

CHECK(Size getNumberOf(const Element* element) const)
	const Element* e = e_ptr->getElement(6);
	Size num3 = e_ptr->getNumberOf(e);
	TEST_EQUAL(num3, 4);
RESULT

CHECK(Size getNumberOfAtoms() const)
	Size num4 = e_ptr->getNumberOfAtoms();
	TEST_EQUAL(num4, 4);
RESULT

CHECK(EmpiricalFormula& operator = (const EmpiricalFormula& rhs))
	EmpiricalFormula ef;
	ef = *e_ptr;
	TEST_EQUAL(*e_ptr == ef, true)
RESULT

CHECK(EmpiricalFormula& operator = (const String& rhs) throw(Exception::ParseError))
	EmpiricalFormula ef;
	ef = "C4";
	TEST_EQUAL(*e_ptr == ef, true)
RESULT


CHECK(EmpiricalFormula& operator += (const EmpiricalFormula& rhs))
	EmpiricalFormula ef("C3");
	ef += ef;
	TEST_EQUAL(ef, "C6")
RESULT

CHECK(EmpiricalFormula& operator += (const String& rhs) throw(Exception::ParseError))
	EmpiricalFormula ef;
	ef += "C";
	TEST_EQUAL(ef, "C")
RESULT

CHECK(EmpiricalFormula operator + (const EmpiricalFormula& rhs) const)
	EmpiricalFormula ef("C2");
	EmpiricalFormula ef2;
	ef2 = ef + ef;
	TEST_EQUAL(ef2, "C4")
RESULT

CHECK(EmpiricalFormula operator + (const String& rhs) const throw(Exception::ParseError))
	EmpiricalFormula ef1("C2");
	EmpiricalFormula ef2;
	ef2 = ef1 + "C2";
	TEST_EQUAL(ef2, "C4")
RESULT

CHECK(EmpiricalFormula& operator -= (const EmpiricalFormula& rhs) throw(Exception::SizeUnderflow))
	EmpiricalFormula ef1("C5H12"), ef2("CH12");
	ef1 -= ef2;
	TEST_EQUAL(*e_ptr == ef1, true)
	TEST_EXCEPTION(Exception::SizeUnderflow, ef1 -= ef2)
RESULT

CHECK(EmpiricalFormula& operator -= (const String& rhs) throw(Exception::ParseError, Exception::SizeUnderflow))
	EmpiricalFormula ef1("C5H12");
	ef1 -= "CH12";
	TEST_EQUAL(*e_ptr == ef1, true)
	TEST_EXCEPTION(Exception::SizeUnderflow, ef1 -= "CH12")
RESULT

CHECK(EmpiricalFormula operator - (const EmpiricalFormula& rhs) const throw(Exception::SizeUnderflow))
	EmpiricalFormula ef1("C5H12"), ef2("CH12");
	EmpiricalFormula ef3, ef4;
	ef3 = ef1 - ef2;
	TEST_EQUAL(*e_ptr == ef3, true)
RESULT

CHECK(EmpiricalFormula operator - (const String& rhs) const throw(Exception::ParseError, Exception::SizeUnderflow))
	EmpiricalFormula ef1("C5H12"), ef2("CH12"), ef4;
	ef4 = ef1 - "CH12";
	TEST_EQUAL(*e_ptr == ef4, true)
	TEST_EXCEPTION(Exception::SizeUnderflow, ef1-"O3")
	TEST_EXCEPTION(Exception::SizeUnderflow, ef1-"C6")
	TEST_EXCEPTION(Exception::SizeUnderflow, ef2-ef1)
	TEST_EXCEPTION(Exception::ParseError, ef1-"BLUBB")
RESULT

CHECK(bool isEmpty() const)
	EmpiricalFormula ef;
	TEST_EQUAL(ef.isEmpty(), true)
	TEST_EQUAL(e_ptr->isEmpty(), false)
RESULT

CHECK(bool hasElement(const String& name) const)
	TEST_EQUAL(e_ptr->hasElement("C"), true)
	TEST_EQUAL(e_ptr->hasElement("N"), false)
RESULT

CHECK(bool hasElement(Size atomic_number) const)
	TEST_EQUAL(e_ptr->hasElement(6), true)
	TEST_EQUAL(e_ptr->hasElement(7), false)
RESULT

CHECK(bool hasElement(const Element* element) const)
	const Element* e = e_ptr->getElement(6);
	TEST_EQUAL(e_ptr->hasElement(e), true)
	e = e_ptr->getElement(1);
	TEST_EQUAL(e_ptr->hasElement(e), false)
RESULT

CHECK(void setCharge(SignedInt charge))
	e_ptr->setCharge(1);
RESULT

CHECK(SignedInt getCharge() const)
	TEST_EQUAL(e_ptr->getCharge(), 1)
RESULT

CHECK(bool isCharged() const)
	TEST_EQUAL(e_ptr->isCharged(), true)
	e_ptr->setCharge(0);
	TEST_EQUAL(e_ptr->isCharged(), false)
RESULT

CHECK(Real getAverageWeight() const)
	EmpiricalFormula ef("C2");
	const Element* e = e_ptr->getElement("C");
	TEST_REAL_EQUAL(ef.getAverageWeight(), e->getAverageWeight() * 2)
RESULT

CHECK(Real getMonoWeight() const)
	EmpiricalFormula ef("C2");
	const Element* e = e_ptr->getElement("C");
	TEST_REAL_EQUAL(ef.getMonoWeight(), e->getMonoWeight() * 2)
RESULT

CHECK(String getString() const)
	EmpiricalFormula ef("C2H5");
	String str = ef.getString();
	TEST_EQUAL(String(str).hasSubstring("H5"), true)
	TEST_EQUAL(String(str).hasSubstring("C2"), true)
RESULT

CHECK(const ElementDB* getElementDB() const)
	const ElementDB* db = 0;
	db = e_ptr->getElementDB();
	TEST_NOT_EQUAL(db, 0);
	TEST_EQUAL(db->getElement("C")->getSymbol(), "C")
RESULT

CHECK(void setElementDB(const String& file_name) throw(Exception::FileNotFound, Exception::ParseError))
	// TODO
RESULT

CHECK((friend std::ostream& operator << (std::ostream&, const EmpiricalFormula&)))
	stringstream ss;
	EmpiricalFormula ef("C2H5");
	ss << ef;
	TEST_EQUAL(String(ss.str()).hasSubstring("H5"), true);
	TEST_EQUAL(String(ss.str()).hasSubstring("C2"), true);
RESULT

CHECK(bool operator != (const EmpiricalFormula& rhs) const)
	EmpiricalFormula ef1("C2H5"), ef2(*e_ptr);
	TEST_EQUAL(ef1 != ef2, true)
	TEST_EQUAL(ef1 != ef1, false)
	ef2.setCharge(1);
	TEST_EQUAL(ef2 != *e_ptr, true)
RESULT

CHECK(bool operator != (const String& rhs) const throw(Exception::ParseError))
	EmpiricalFormula ef1("C2H5");
	TEST_EQUAL(ef1 != "C2", true)
	TEST_EQUAL(ef1 != "C2H5", false)
RESULT

CHECK(bool operator == (const EmpiricalFormula& rhs) const)
	EmpiricalFormula ef1("C2H5"), ef2(*e_ptr);
	TEST_EQUAL(ef1 == ef2, false)
	TEST_EQUAL(ef1 == ef1, true)
	ef2.setCharge(1);
	TEST_EQUAL(ef2 == *e_ptr, false)
RESULT

CHECK(bool operator == (const String& rhs) const throw(Exception::ParseError))
	EmpiricalFormula ef1("C2H5");
	TEST_EQUAL(ef1 == "C2", false)
	TEST_EQUAL(ef1 == "C2H5", true)
RESULT

CHECK(ConstIterator begin() const)
 // TODO
RESULT

CHECK(ConstIterator end() const)
 // TODO
RESULT

CHECK(IsotopeDistribution getIsotopeDistribution(Size max_depth = 20) const)
	// TODO
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
