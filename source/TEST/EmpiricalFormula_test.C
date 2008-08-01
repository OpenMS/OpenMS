// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

CHECK(EmpiricalFormula(const String& rhs))
	e_ptr = new EmpiricalFormula("C4");
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(EmpiricalFormula(const EmpiricalFormula& rhs))
	EmpiricalFormula ef(*e_ptr);
	TEST_EQUAL(ef == *e_ptr, true)
RESULT

CHECK((EmpiricalFormula(UInt number, const Element* element, Int charge=0)))
	EmpiricalFormula ef(4, e_ptr->getElement("C"));
	TEST_EQUAL(ef == *e_ptr, true)
	TEST_EQUAL(ef.getCharge(), 0)
RESULT

CHECK(const Element* getElement(UInt atomic_number) const)
	const Element* e = e_ptr->getElement(6);
	TEST_EQUAL(e->getSymbol(), "C")
RESULT

CHECK(const Element* getElement(const String& name) const)
	const Element* e = e_ptr->getElement("C");
	TEST_EQUAL(e->getSymbol(), "C")
RESULT

CHECK(UInt getNumberOf(UInt atomic_number) const)
	UInt num1 = e_ptr->getNumberOf(6);
	TEST_EQUAL(num1, 4);
RESULT

CHECK(UInt getNumberOf(const String& name) const)
	UInt num2 = e_ptr->getNumberOf("C");
	TEST_EQUAL(num2, 4);
RESULT

CHECK(UInt getNumberOf(const Element* element) const)
	const Element* e = e_ptr->getElement(6);
	UInt num3 = e_ptr->getNumberOf(e);
	TEST_EQUAL(num3, 4);
RESULT

CHECK(UInt getNumberOfAtoms() const)
	UInt num4 = e_ptr->getNumberOfAtoms();
	TEST_EQUAL(num4, 4);
RESULT

CHECK(EmpiricalFormula& operator = (const EmpiricalFormula& rhs))
	EmpiricalFormula ef;
	ef = *e_ptr;
	TEST_EQUAL(*e_ptr == ef, true)
RESULT

CHECK(EmpiricalFormula& operator = (const String& rhs))
	EmpiricalFormula ef;
	ef = "C4";
	TEST_EQUAL(*e_ptr == ef, true)
	TEST_EXCEPTION(Exception::ParseError, ef = "2C4")
RESULT


CHECK(EmpiricalFormula& operator += (const EmpiricalFormula& rhs))
	EmpiricalFormula ef("C3");
	ef += ef;
	TEST_EQUAL(ef, "C6")
RESULT

CHECK(EmpiricalFormula& operator += (const String& rhs))
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

CHECK(EmpiricalFormula operator + (const String& rhs) const)
	EmpiricalFormula ef1("C2");
	EmpiricalFormula ef2;
	ef2 = ef1 + "C2";
	TEST_EQUAL(ef2, "C4")
RESULT

CHECK(EmpiricalFormula& operator -= (const EmpiricalFormula& rhs))
	EmpiricalFormula ef1("C5H12"), ef2("CH12");
	ef1 -= ef2;
	TEST_EQUAL(*e_ptr == ef1, true)
	TEST_EXCEPTION(Exception::SizeUnderflow, ef1 -= ef2)
RESULT

CHECK(EmpiricalFormula& operator -= (const String& rhs))
	EmpiricalFormula ef1("C5H12");
	ef1 -= "CH12";
	TEST_EQUAL(*e_ptr == ef1, true)
	TEST_EXCEPTION(Exception::SizeUnderflow, ef1 -= "CH12")
RESULT

CHECK(EmpiricalFormula operator - (const EmpiricalFormula& rhs) const)
	EmpiricalFormula ef1("C5H12"), ef2("CH12");
	EmpiricalFormula ef3, ef4;
	ef3 = ef1 - ef2;
	TEST_EQUAL(*e_ptr == ef3, true)
RESULT

CHECK(EmpiricalFormula operator - (const String& rhs) const)
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

CHECK(bool hasElement(UInt atomic_number) const)
	TEST_EQUAL(e_ptr->hasElement(6), true)
	TEST_EQUAL(e_ptr->hasElement(7), false)
RESULT

CHECK(bool hasElement(const Element* element) const)
	const Element* e = e_ptr->getElement(6);
	TEST_EQUAL(e_ptr->hasElement(e), true)
	e = e_ptr->getElement(1);
	TEST_EQUAL(e_ptr->hasElement(e), false)
RESULT

CHECK(void setCharge(Int charge))
	e_ptr->setCharge(1);
	NOT_TESTABLE // will be tested in next check
RESULT

CHECK(Int getCharge() const)
	TEST_EQUAL(e_ptr->getCharge(), 1)
RESULT

CHECK(bool isCharged() const)
	TEST_EQUAL(e_ptr->isCharged(), true)
	e_ptr->setCharge(0);
	TEST_EQUAL(e_ptr->isCharged(), false)
RESULT

CHECK(DoubleReal getAverageWeight() const)
	EmpiricalFormula ef("C2");
	const Element* e = e_ptr->getElement("C");
	TEST_REAL_EQUAL(ef.getAverageWeight(), e->getAverageWeight() * 2)
RESULT

CHECK(DoubleReal getMonoWeight() const)
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

CHECK([EXTRA](friend std::ostream& operator << (std::ostream&, const EmpiricalFormula&)))
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

CHECK(bool operator != (const String& rhs) const)
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

CHECK(bool operator == (const String& rhs) const)
	EmpiricalFormula ef1("C2H5");
	TEST_EQUAL(ef1 == "C2", false)
	TEST_EQUAL(ef1 == "C2H5", true)
RESULT

CHECK(ConstIterator begin() const)
	EmpiricalFormula ef("C6H12O6");
	Map<String, UInt> formula;
	formula["C"] = 6;
	formula["H"] = 12;
	formula["O"] = 6;
	for (EmpiricalFormula::ConstIterator it = ef.begin(); it != ef.end(); ++it)
	{
		TEST_EQUAL(it->second, formula[it->first->getSymbol()])
	}

RESULT

CHECK(ConstIterator end() const)
	NOT_TESTABLE
RESULT

CHECK(IsotopeDistribution getIsotopeDistribution(UInt max_depth) const)
	EmpiricalFormula ef("C");
	IsotopeDistribution iso = ef.getIsotopeDistribution(20);
	double result[] = { 0.9893, 0.0107};
	UInt i = 0;
	for (IsotopeDistribution::ConstIterator it = iso.begin(); it != iso.end(); ++it, ++i)
	{
		TEST_REAL_EQUAL(it->second, result[i])
	}
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

