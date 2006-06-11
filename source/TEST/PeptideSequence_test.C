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
// $Id: PeptideSequence_test.C,v 1.3 2006/03/28 12:53:13 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/PeptideSequence.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ResidueDB, "$Id: PeptideSequence_test.C,v 1.3 2006/03/28 12:53:13 marc_sturm Exp $")

/////////////////////////////////////////////////////////////

PeptideSequence* p_ptr = 0;
CHECK(PeptideSequence())
	p_ptr = new PeptideSequence;
	TEST_NOT_EQUAL(p_ptr, 0)
RESULT

CHECK(~PeptideSequence())
	delete p_ptr;
RESULT

CHECK(operator +)
	PeptideSequence seq1("AAACCC");
	PeptideSequence seq2("AAADDD");
	TEST_EQUAL(seq1+seq2 == "AAACCCAAADDD", true)
	TEST_EQUAL(seq1+seq2 != "AAACCCAAADDD", false)
	TEST_EQUAL(seq1+seq2 == "AAA", false)
	PeptideSequence seq3("AlaAlaAlaAspAspAsp");
	TEST_EQUAL(seq2 == seq3, true)
	TEST_EQUAL(seq1.hasPrefix("Ala"), true)
	TEST_EQUAL(seq1.hasSubsequence("Ala"), true)
	TEST_EXCEPTION(Exception::ParseError, PeptideSequence("BZ"))
	TEST_EXCEPTION(Exception::ParseError, PeptideSequence("M(Ox"))
	TEST_EXCEPTION(Exception::ParseError, PeptideSequence("A(Ox)"))
	PeptideSequence seq4("M(Ox)MM");
	TEST_EQUAL(seq4.hasPrefix(PeptideSequence("M(Ox)")), true)
	
	TEST_EQUAL(seq1.getSuffix(3) == "CCC", true)
	TEST_EQUAL(seq1.getPrefix(3) == "AAA", true)
	TEST_EQUAL(seq1.getPrefix(1) == "A", true)
	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
