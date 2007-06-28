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

CHECK((const Residue* getResidue(Int index) const throw(Exception::IndexUnderflow, Exception::IndexOverflow)))
	AASequence seq(String("ACDEF"));
	Int sint(2);
	TEST_EQUAL(seq.getResidue(sint)->getOneLetterCode(), "D")
	TEST_EXCEPTION(Exception::IndexUnderflow, seq.getResidue(-3))
	TEST_EXCEPTION(Exception::IndexOverflow, seq.getResidue(1000))
RESULT

CHECK(const Residue* getResidue(UInt index) const throw(Exception::IndexOverflow))
	AASequence seq("ACDEF");
	UInt unsignedint(2);
	TEST_EQUAL(seq.getResidue(unsignedint)->getOneLetterCode(), "D")
	TEST_EXCEPTION(Exception::IndexOverflow, seq.getResidue(1000))
RESULT

CHECK((EmpiricalFormula getFormula(Residue::ResidueType type = Residue::Full, Int charge=0) const))
	AASequence seq("ACDEF");
	TEST_EQUAL(seq.getFormula(), "SC26O12N5H35")
	TEST_EQUAL(seq.getFormula(Residue::Full, 1), "SC26O12N5H36")
	TEST_EQUAL(seq.getFormula(Residue::BIon, 0), "SC26O11N5H33")
RESULT

CHECK((DoubleReal getAverageWeight(Residue::ResidueType type = Residue::Full, Int charge=0) const))
	AASequence seq("DFPIANGER");
	PRECISION(0.01)
	TEST_REAL_EQUAL(seq.getAverageWeight(), double(1018.08))
	TEST_REAL_EQUAL(seq.getAverageWeight(Residue::YIon, 1), double(1019.09))
RESULT

CHECK((DoubleReal getMonoWeight(Residue::ResidueType type = Residue::Full, Int charge=0) const))
  AASequence seq("DFPIANGER");
	PRECISION(0.01)
	TEST_REAL_EQUAL(seq.getMonoWeight(), double(1017.49))
	TEST_REAL_EQUAL(seq.getMonoWeight(Residue::YIon, 1), double(1018.5))
RESULT

CHECK((HashMap<const EmpiricalFormula*, UInt> getNeutralLosses() const))
	AASequence seq("DFPIANGER");
  HashMap<const EmpiricalFormula*, UInt> losses = seq.getNeutralLosses();
	TEST_EQUAL(losses.size(), 10)
RESULT

CHECK(const Residue* operator [] (Int index) const throw(Exception::IndexUnderflow, Exception::IndexOverflow))
  AASequence seq("DFPIANGER");
	Int index = 0;
	TEST_EQUAL(seq[index]->getOneLetterCode(), "D")
	index = -1;
	TEST_EXCEPTION(Exception::IndexUnderflow, seq[index])
	index = 20;
	TEST_EXCEPTION(Exception::IndexOverflow, seq[index])
RESULT

CHECK(const Residue* operator [] (UInt index) const throw(Exception::IndexOverflow))
	AASequence seq("DFPIANGER");
  UInt index = 0;
  TEST_EQUAL(seq[index]->getOneLetterCode(), "D")
  index = 20;
  TEST_EXCEPTION(Exception::IndexOverflow, seq[index])
RESULT

CHECK(AASequence operator + (const AASequence& peptide) const)
  AASequence seq1("DFPIANGER"), seq2("DFP"), seq3("IANGER");
	TEST_EQUAL(seq1, seq2 + seq3);
RESULT

CHECK(AASequence operator + (const String& peptide) const throw(Exception::ParseError))
  AASequence seq1("DFPIANGER"), seq2("DFP"); 
	String seq3("IANGER"), seq4("BLUBB");
	TEST_EQUAL(seq1, seq2 + seq3)
	TEST_EXCEPTION(Exception::ParseError, seq2 + seq4)
RESULT

CHECK(AASequence& operator += (const AASequence&))
  AASequence seq1("DFPIANGER"), seq2("DFP"), seq3("IANGER");
	seq2 += seq3;
	TEST_EQUAL(seq1, seq2)
RESULT

CHECK(AASequence& operator += (const String&) throw(Exception::ParseError))
  AASequence seq1("DFPIANGER"), seq2("DFP");
	String seq3("IANGER"), seq4("BLUBB");
	seq2 += seq3;
	TEST_EQUAL(seq1, seq2)
	TEST_EXCEPTION(Exception::ParseError, seq2 += seq4)
RESULT

CHECK((void setResidueDB(ResidueDB* res_db=0)))
	ResidueDB* res_db = new ResidueDB();
	AASequence seq;
	seq.setResidueDB(res_db);
RESULT

CHECK(UInt size() const)
  AASequence seq1("DFPIANGER");
	TEST_EQUAL(seq1.size(), 9)
RESULT

CHECK(AASequence getPrefix(UInt index) const throw(Exception::IndexOverflow))
  AASequence seq1("DFPIANGER"), seq2("DFP"), seq3("DFPIANGER");
	TEST_EQUAL(seq2, seq1.getPrefix(3));
	TEST_EQUAL(seq3, seq1.getPrefix(9));
	TEST_EXCEPTION(Exception::IndexOverflow, seq1.getPrefix(10))
RESULT

CHECK(AASequence getSuffix(UInt index) const throw(Exception::IndexOverflow))
  AASequence seq1("DFPIANGER"), seq2("GER"), seq3("DFPIANGER");
	TEST_EQUAL(seq2, seq1.getSuffix(3));
	TEST_EQUAL(seq3, seq1.getSuffix(9));
	TEST_EXCEPTION(Exception::IndexOverflow, seq1.getSuffix(10))
RESULT

CHECK(AASequence getSubsequence(UInt index, UInt number) const throw(Exception::IndexOverflow))
  AASequence seq1("DFPIANGER"), seq2("IAN"), seq3("DFPIANGER");
	TEST_EQUAL(seq2, seq1.getSubsequence(3, 3))
	TEST_EQUAL(seq3, seq1.getSubsequence(0, 9))
	TEST_EXCEPTION(Exception::IndexOverflow, seq1.getSubsequence(0, 10))
RESULT

CHECK(bool has(const Residue* residue) const)
  AASequence seq("DFPIANGER");
	TEST_EQUAL(seq.has(seq[0]), true)
	Residue* res = 0;
	TEST_NOT_EQUAL(seq.has(res), true)
RESULT

CHECK(bool has(const String& name) const)
  AASequence seq("DFPIANGER");
	TEST_EQUAL(seq.has("D"), true)
	TEST_EQUAL(seq.has("N"), true)
	TEST_NOT_EQUAL(seq.has("Q"), true)
RESULT

CHECK(bool hasSubsequence(const AASequence& peptide) const)
  AASequence seq1("DFPIANGER"), seq2("IANG"), seq3("AIN");
	TEST_EQUAL(seq1.hasSubsequence(seq2), true)
	TEST_EQUAL(seq1.hasSubsequence(seq3), false)
RESULT

CHECK(bool hasSubsequence(const String& peptide) const throw(Exception::ParseError))
  AASequence seq1("DFPIANGER");
	String seq2("IANG"), seq3("AIN");
	TEST_EQUAL(seq1.hasSubsequence(seq2), true)
	TEST_EQUAL(seq1.hasSubsequence(seq3), false)
RESULT

CHECK(bool hasPrefix(const AASequence& peptide) const)
  AASequence seq1("DFPIANGER"), seq2("DFP"), seq3("AIN");
	TEST_EQUAL(seq1.hasPrefix(seq2), true)
	TEST_EQUAL(seq1.hasPrefix(seq3), false)
RESULT

CHECK(bool hasPrefix(const String& peptide) const throw(Exception::ParseError))
  AASequence seq1("DFPIANGER");
	String seq2("DFP"), seq3("AIN"), seq4("BLUBB");
	TEST_EQUAL(seq1.hasPrefix(seq2), true)
	TEST_EQUAL(seq1.hasPrefix(seq3), false)
TEST_EXCEPTION(Exception::ParseError, seq1.hasPrefix(seq4))
RESULT

CHECK(bool hasSuffix(const AASequence& peptide) const)
  AASequence seq1("DFPIANGER"), seq2("GER"), seq3("AIN");
  TEST_EQUAL(seq1.hasSuffix(seq2), true)
  TEST_EQUAL(seq1.hasSuffix(seq3), false) 
RESULT

CHECK(bool hasSuffix(const String& peptide) const throw(Exception::ParseError))
  AASequence seq1("DFPIANGER");
  String seq2("GER"), seq3("AIN"), seq4("BLUBB");
  TEST_EQUAL(seq1.hasSuffix(seq2), true)
  TEST_EQUAL(seq1.hasSuffix(seq3), false)
	TEST_EXCEPTION(Exception::ParseError, seq1.hasSuffix(seq4))
RESULT

CHECK(ConstIterator begin() const)
	String result[] = { "D", "F", "P", "I", "A", "N", "G", "E", "R" };
	AASequence seq("DFPIANGER");
	UInt i = 0;
	for (AASequence::ConstIterator it = seq.begin(); it != seq.end(); ++it, ++i)
	{
		TEST_EQUAL((*it)->getOneLetterCode(), result[i])
	}
RESULT

CHECK(ConstIterator end() const)
	// testet above
RESULT

CHECK(Iterator begin())
	String result[] = { "D", "F", "P", "I", "A", "N", "G", "E", "R" };
  AASequence seq("DFPIANGER");
  UInt i = 0;
  for (AASequence::ConstIterator it = seq.begin(); it != seq.end(); ++it, ++i)
  {
    TEST_EQUAL((*it)->getOneLetterCode(), result[i])
  }
RESULT
		  
CHECK(Iterator end())
	// tested above
RESULT

//CHECK(friend std::ostream& operator << (std::ostream& os, const AASequence& peptide))
//  // TODO
//RESULT

//CHECK(friend std::istream& operator >> (std::istream& is, const AASequence& peptide))
//  // TODO
//RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
