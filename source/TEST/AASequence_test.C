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
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ResidueDB, "$Id$")

/////////////////////////////////////////////////////////////

AASequence* ptr = 0;
CHECK(AASequence())
	ptr = new AASequence();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~AASequence())
	delete ptr;
RESULT

CHECK(AASequence(const AASequence& rhs))
	AASequence seq;
	seq = String("AAA");
	AASequence seq2(seq);
	TEST_EQUAL(seq, seq2)
RESULT

CHECK(AASequence(const String& rhs) throw(Exception::ParseError))
	AASequence seq;
	seq = String("AAA");
	AASequence seq2("AAA");
	TEST_EQUAL(seq, seq2);
RESULT

CHECK(AASequence& operator = (const AASequence& rhs))
	AASequence seq("AAA");
	AASequence seq2;
	seq2 = String("AAA");
	TEST_EQUAL(seq, seq2)
RESULT

CHECK(AASequence(ResidueDB* res_db))
	ResidueDB* res_db = new ResidueDB();
	AASequence seq(res_db);
RESULT

CHECK(AASequence(ConstIterator begin, ConstIterator end))
	AASequence seq("ACDEFGHIKLMN");
	AASequence seq2(seq.begin(), seq.end() - 4);
	AASequence seq3("ACDEFGHI");
	TEST_EQUAL(seq2, seq3);
RESULT

CHECK(bool operator == (const AASequence&) const)
	AASequence seq1("ACDEF");
	AASequence seq2("ADCEF");
	TEST_EQUAL(seq1 == seq1, true)
	TEST_EQUAL(seq1 == seq2, false)
RESULT

CHECK(bool operator == (const String&) const throw(Exception::ParseError))
	AASequence seq("ACDEF");
	TEST_EQUAL(seq == "ACDEF", true)
	TEST_EQUAL(seq == "ADCEF", false)
RESULT

CHECK(bool operator != (const AASequence&) const)
  AASequence seq1("ACDEF");
  AASequence seq2("ADCEF");
  TEST_EQUAL(seq1 != seq1, false)
  TEST_EQUAL(seq1 != seq2, true)
RESULT

CHECK(bool operator != (const String&) const throw(Exception::ParseError))
  AASequence seq("ACDEF");
  TEST_EQUAL(seq != "ACDEF", false)
  TEST_EQUAL(seq != "ADCEF", true)
RESULT

CHECK((const Residue* getResidue(SignedInt index) const throw(Exception::IndexUnderflow, Exception::IndexOverflow)))
	AASequence seq(String("ACDEF"));
	SignedInt sint(2);
	TEST_EQUAL(seq.getResidue(sint)->getOneLetterCode(), "D")
	TEST_EXCEPTION(Exception::IndexUnderflow, seq.getResidue(-3))
	TEST_EXCEPTION(Exception::IndexOverflow, seq.getResidue(1000))
RESULT

CHECK(const Residue* getResidue(UnsignedInt index) const throw(Exception::IndexOverflow))
	AASequence seq("ACDEF");
	UnsignedInt unsignedint(2);
	TEST_EQUAL(seq.getResidue(unsignedint)->getOneLetterCode(), "D")
	TEST_EXCEPTION(Exception::IndexOverflow, seq.getResidue(1000))
RESULT

CHECK(EmpiricalFormula getFormula(Residue::ResidueType type = Residue::Full, SignedInt charge = 0) const)
	AASequence seq("ACDEF");
	TEST_EQUAL(seq.getFormula(), "SC26O12N5H35")
	TEST_EQUAL(seq.getFormula(Residue::Full, 1), "SC26O12N5H36")
	TEST_EQUAL(seq.getFormula(Residue::BIon, 0), "SC26O11N5H33")
RESULT

CHECK(Real getAverageWeight(Residue::ResidueType type = Residue::Full, SignedInt charge = 0) const)

RESULT

CHECK(Real getMonoWeight(Residue::ResidueType type = Residue::Full, SignedInt charge = 0) const)

RESULT

CHECK((HashMap<const EmpiricalFormula*, Size> getNeutralLosses() const))

RESULT

CHECK(const Residue* operator [] (SignedInt index) const throw(Exception::IndexUnderflow, Exception::IndexOverflow))

RESULT

CHECK(const Residue* operator [] (UnsignedInt index) const throw(Exception::IndexOverflow))

RESULT

CHECK(AASequence operator + (const AASequence& peptide) const)

RESULT

CHECK(AASequence operator + (const String& peptide) const throw(Exception::ParseError))

RESULT

CHECK(AASequence& operator += (const AASequence&))

RESULT

CHECK(AASequence& operator += (const String&) throw(Exception::ParseError))

RESULT

CHECK(void setResidueDB(ResidueDB* res_db = 0))
	// TODO
RESULT

CHECK(Size size() const)

RESULT

CHECK(AASequence getPrefix(Size index) const throw(Exception::IndexOverflow))

RESULT

CHECK(AASequence getSuffix(Size index) const throw(Exception::IndexOverflow))

RESULT

CHECK(AASequence getSubsequence(Size index, Size number) const throw(Exception::IndexOverflow))

RESULT

CHECK(bool has(const Residue* residue) const)

RESULT

CHECK(bool has(const String& name) const)

RESULT

CHECK(bool hasSubsequence(const AASequence& peptide) const)

RESULT

CHECK(bool hasSubsequence(const String& peptide) const throw(Exception::ParseError))

RESULT

CHECK(bool hasPrefix(const AASequence& peptide) const)

RESULT

CHECK(bool hasPrefix(const String& peptide) const throw(Exception::ParseError))

RESULT

CHECK(bool hasSuffix(const AASequence& peptide) const)

RESULT

CHECK(bool hasSuffix(const String& peptide) const throw(Exception::ParseError))

RESULT

CHECK(ConstIterator begin() const)
	// TODO
RESULT

CHECK(ConstIterator end() const)
	// TODO
RESULT

CHECK(Iterator begin())
	// TODO
RESULT
		  
CHECK(Iterator end())
	// TODO
RESULT

CHECK(friend std::ostream& operator << (std::ostream& os, const AASequence& peptide))

RESULT

CHECK(friend std::istream& operator >> (std::istream& is, const AASequence& peptide))

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
