// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: EmpiricalFormula_test.C,v 1.3 2006/03/28 12:53:13 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ElementDB, "$Id: EmpiricalFormula_test.C,v 1.3 2006/03/28 12:53:13 marc_sturm Exp $")

/////////////////////////////////////////////////////////////

EmpiricalFormula* e_ptr = 0;
CHECK(EmpiricalFormula())
	e_ptr = new EmpiricalFormula;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~EmpiricalFormula())
	delete e_ptr;
RESULT

CHECK(EmpiricalFormula(const String&))
	e_ptr = new EmpiricalFormula("C4");
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(EmpiricalFormula(const EmpiricalFormula&))
	EmpiricalFormula ef(*e_ptr);
	TEST_EQUAL(ef == *e_ptr, true)
RESULT

CHECK(EmpiricalFormula(Size number, const Element*))
	EmpiricalFormula ef(4, e_ptr->getElement("C"));
	TEST_EQUAL(ef == *e_ptr, true)
RESULT

CHECK(getElement())
	const Element * e1 = e_ptr->getElement(6);
	const Element * e2 = e_ptr->getElement("C");
	TEST_EQUAL(e1, e2)
	TEST_NOT_EQUAL(e1, 0)
RESULT

CHECK(getNumberOf())
	const Element * e = e_ptr->getElement(6);
	Size num1 = e_ptr->getNumberOf(6);
	TEST_EQUAL(num1, 4);
	Size num2 = e_ptr->getNumberOf("C");
	TEST_EQUAL(num2, 4);
	Size num3 = e_ptr->getNumberOf(e);
	TEST_EQUAL(num3, 4);
	Size num4 = e_ptr->getNumberOfAtoms();
	TEST_EQUAL(num4, 4);
RESULT

CHECK(operator =)
	EmpiricalFormula ef1, ef2;
	ef1 = *e_ptr;
	ef2 = "C4";
	TEST_EQUAL(*e_ptr == ef1, true)
	TEST_EQUAL(*e_ptr == ef2, true)
RESULT

CHECK(operator +=)
	EmpiricalFormula ef1("C3");
	EmpiricalFormula ef2;
	ef2 += ef1;
	ef2 += "C";
	TEST_EQUAL(*e_ptr == ef2, true)
	TEST_EQUAL(*e_ptr != ef1, true)
RESULT

CHECK(operator +)
	EmpiricalFormula ef1("C2");
	EmpiricalFormula ef2, ef3;
	ef2 = ef1 + ef1;
	ef3 = ef1 + "C2";
	TEST_EQUAL(*e_ptr == ef2, true)
	TEST_EQUAL(*e_ptr == ef3, true)
RESULT

CHECK(operator -=)
	EmpiricalFormula ef1("C5H12"), ef2("CH12");
	ef1 -= ef2;
	TEST_EQUAL(*e_ptr == ef1, true)
	TEST_EXCEPTION(Exception::SizeUnderflow, ef1 -= ef2)

	ef1 = "C5H12";
	ef1 -= "CH12";
	TEST_EQUAL(*e_ptr == ef1, true)
	TEST_EXCEPTION(Exception::SizeUnderflow, ef1 -= "CH12")
RESULT

CHECK(operator -)
	EmpiricalFormula ef1("C5H12"), ef2("CH12");
	EmpiricalFormula ef3, ef4;
	ef3 = ef1 - ef2;
	ef4 = ef1 - "CH12";
	TEST_EQUAL(*e_ptr == ef3, true)
	TEST_EQUAL(*e_ptr == ef4, true)
	TEST_EXCEPTION(Exception::SizeUnderflow, ef1-"O3")
	TEST_EXCEPTION(Exception::SizeUnderflow, ef1-"C6")
	TEST_EXCEPTION(Exception::SizeUnderflow, ef2-ef1)
	TEST_EXCEPTION(Exception::ParseError, ef1-"BLUBB")
RESULT

CHECK(isEmpty())
	EmpiricalFormula ef;
	TEST_EQUAL(ef.isEmpty(), true)
	TEST_EQUAL(e_ptr->isEmpty(), false)
RESULT

CHECK(hasElement())
	TEST_EQUAL(e_ptr->hasElement("C"), true)
	TEST_EQUAL(e_ptr->hasElement("N"), false)
	TEST_EQUAL(e_ptr->hasElement(6), true)
	TEST_EQUAL(e_ptr->hasElement(7), false)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
