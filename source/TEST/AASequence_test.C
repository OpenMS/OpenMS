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

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ResidueDB, "$Id$")

/////////////////////////////////////////////////////////////

AASequence* p_ptr = 0;
CHECK(AASequence())
	p_ptr = new AASequence;
	TEST_NOT_EQUAL(p_ptr, 0)
RESULT

CHECK(~AASequence())
	delete p_ptr;
RESULT

CHECK(operator +)
	AASequence seq1("AAACCC");
	AASequence seq2("AAADDD");
	TEST_EQUAL(seq1+seq2 == "AAACCCAAADDD", true)
	TEST_EQUAL(seq1+seq2 != "AAACCCAAADDD", false)
	TEST_EQUAL(seq1+seq2 == "AAA", false)
	AASequence seq3("AlaAlaAlaAspAspAsp");
	TEST_EQUAL(seq2 == seq3, true)
	TEST_EQUAL(seq1.hasPrefix("Ala"), true)
	TEST_EQUAL(seq1.hasSubsequence("Ala"), true)
	TEST_EXCEPTION(Exception::ParseError, AASequence("BZ"))
	TEST_EXCEPTION(Exception::ParseError, AASequence("M(Ox"))
	TEST_EXCEPTION(Exception::ParseError, AASequence("A(Ox)"))
	AASequence seq4("M(Ox)MM");
	TEST_EQUAL(seq4.hasPrefix(AASequence("M(Ox)")), true)
	
	TEST_EQUAL(seq1.getSuffix(3) == "CCC", true)
	TEST_EQUAL(seq1.getPrefix(3) == "AAA", true)
	TEST_EQUAL(seq1.getPrefix(1) == "A", true)
	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
