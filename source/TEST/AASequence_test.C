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

CHECK(AASequence(const String& rhs))
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

CHECK(AASequence(ConstIterator begin, ConstIterator end))
	const AASequence seq("ACDEFGHIKLMN");
	AASequence seq2(seq.begin(), seq.end() - 4);
	AASequence seq3("ACDEFGHI");
	TEST_EQUAL(seq2, seq3);
RESULT

CHECK(bool operator == (const char* rhs) const)
  AASequence seq1("(MOD:00051)DFPIANGER");
  AASequence seq2("DFPIANGER");
  TEST_EQUAL(seq2 == "DFPIANGER", true)
  TEST_EQUAL(seq1 == "(MOD:00051)DFPIANGER", true)

  AASequence seq3("DFPIANGER(MOD:00177)");
  AASequence seq4("DFPIANGER(ArgN)");
  TEST_EQUAL(seq3 == "DFPIANGER", false)
  TEST_EQUAL(seq3 == "DFPIANGER(MOD:00177)", true)
  TEST_EQUAL(seq4 == "DFPIANGER(ArgN)", true)
  TEST_EQUAL(seq4 == "DFPIANGER", false)

  AASequence seq5("DFBIANGER");
  TEST_EQUAL(seq5 == "DFPIANGER", false)
  TEST_EQUAL(seq5 == "DFBIANGER", true)	
RESULT

CHECK(bool operator == (const String& rhs) const)
  AASequence seq1("(MOD:00051)DFPIANGER");
  AASequence seq2("DFPIANGER");
  TEST_EQUAL(seq2 == String("DFPIANGER"), true)
  TEST_EQUAL(seq1 == String("(MOD:00051)DFPIANGER"), true)

  AASequence seq3("DFPIANGER(MOD:00177)");
  AASequence seq4("DFPIANGER(ArgN)");
  TEST_EQUAL(seq3 == String("DFPIANGER"), false)
  TEST_EQUAL(seq3 == String("DFPIANGER(MOD:00177)"), true)
  TEST_EQUAL(seq4 == String("DFPIANGER(ArgN)"), true)
  TEST_EQUAL(seq4 == String("DFPIANGER"), false)

  AASequence seq5("DFBIANGER");
  TEST_EQUAL(seq5 == String("DFPIANGER"), false)
  TEST_EQUAL(seq5 == String("DFBIANGER"), true)
RESULT

CHECK(bool operator == (const AASequence& rhs) const)
  AASequence seq1("(MOD:00051)DFPIANGER");
  AASequence seq2("DFPIANGER");
  TEST_EQUAL(seq2 == AASequence("DFPIANGER"), true)
  TEST_EQUAL(seq1 == AASequence("(MOD:00051)DFPIANGER"), true)

  AASequence seq3("DFPIANGER(MOD:00177)");
  AASequence seq4("DFPIANGER(ArgN)");
  TEST_EQUAL(seq3 == AASequence("DFPIANGER"), false)
  TEST_EQUAL(seq3 == AASequence("DFPIANGER(MOD:00177)"), true)
  TEST_EQUAL(seq4 == AASequence("DFPIANGER(ArgN)"), true)
  TEST_EQUAL(seq4 == AASequence("DFPIANGER"), false)

  AASequence seq5("DFBIANGER");
  TEST_EQUAL(seq5 == AASequence("DFPIANGER"), false)
  TEST_EQUAL(seq5 == AASequence("DFBIANGER"), true)
RESULT

CHECK((const Residue& getResidue(Int index) const))
	AASequence seq(String("ACDEF"));
	Int sint(2);
	TEST_EQUAL(seq.getResidue(sint).getOneLetterCode(), "D")
	TEST_EXCEPTION(Exception::IndexUnderflow, seq.getResidue(-3))
	TEST_EXCEPTION(Exception::IndexOverflow, seq.getResidue(1000))
RESULT

CHECK(const Residue& getResidue(UInt index) const)
	AASequence seq("ACDEF");
	UInt unsignedint(2);
	TEST_EQUAL(seq.getResidue(unsignedint).getOneLetterCode(), "D")
	TEST_EXCEPTION(Exception::IndexOverflow, seq.getResidue(1000))
RESULT

CHECK((EmpiricalFormula getFormula(Residue::ResidueType type = Residue::Full, Int charge=0) const))
	AASequence seq("ACDEF");
	TEST_EQUAL(seq.getFormula(), EmpiricalFormula("O10SH33N5C24"))
	TEST_EQUAL(seq.getFormula(Residue::Full, 1), EmpiricalFormula("O10SH34N5C24"))
	TEST_EQUAL(seq.getFormula(Residue::BIon, 0), EmpiricalFormula("O9SH31N5C24"))
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

CHECK((Map<const EmpiricalFormula*, UInt> getNeutralLosses() const))
	AASequence seq("DFPIANGER");
  Map<const EmpiricalFormula*, UInt> losses = seq.getNeutralLosses();
	TEST_EQUAL(losses.size(), 10)
RESULT

CHECK(const Residue& operator [] (Int index) const)
  AASequence seq("DFPIANGER");
	Int index = 0;
	TEST_EQUAL(seq[index].getOneLetterCode(), "D")
	index = -1;
	TEST_EXCEPTION(Exception::IndexUnderflow, seq[index])
	index = 20;
	TEST_EXCEPTION(Exception::IndexOverflow, seq[index])
RESULT

CHECK(const Residue& operator [] (UInt index) const)
	AASequence seq("DFPIANGER");
  UInt index = 0;
  TEST_EQUAL(seq[index].getOneLetterCode(), "D")
  index = 20;
  TEST_EXCEPTION(Exception::IndexOverflow, seq[index])
RESULT

CHECK(AASequence operator + (const AASequence& peptide) const)
  AASequence seq1("DFPIANGER"), seq2("DFP"), seq3("IANGER");
	TEST_EQUAL(seq1, seq2 + seq3);
RESULT

CHECK(AASequence operator + (const String& peptide) const)
  AASequence seq1("DFPIANGER"), seq2("DFP"); 
	String seq3("IANGER"), seq4("BLUBB");
	TEST_EQUAL(seq1, seq2 + seq3)
RESULT

CHECK(AASequence& operator += (const AASequence&))
  AASequence seq1("DFPIANGER"), seq2("DFP"), seq3("IANGER");
	seq2 += seq3;
	TEST_EQUAL(seq1, seq2)
RESULT

CHECK(AASequence& operator += (const String&))
  AASequence seq1("DFPIANGER"), seq2("DFP");
	String seq3("IANGER"), seq4("BLUBB");
	seq2 += seq3;
	TEST_EQUAL(seq1, seq2)
RESULT

CHECK(UInt size() const)
  AASequence seq1("DFPIANGER");
	TEST_EQUAL(seq1.size(), 9)
RESULT

CHECK(AASequence getPrefix(UInt index) const)
  AASequence seq1("DFPIANGER"), seq2("DFP"), seq3("DFPIANGER");
	TEST_EQUAL(seq2, seq1.getPrefix(3));
	TEST_EQUAL(seq3, seq1.getPrefix(9));
	TEST_EXCEPTION(Exception::IndexOverflow, seq1.getPrefix(10))
RESULT

CHECK(AASequence getSuffix(UInt index) const)
  AASequence seq1("DFPIANGER"), seq2("GER"), seq3("DFPIANGER");
	TEST_EQUAL(seq2, seq1.getSuffix(3));
	TEST_EQUAL(seq3, seq1.getSuffix(9));
	TEST_EXCEPTION(Exception::IndexOverflow, seq1.getSuffix(10))
RESULT

CHECK(AASequence getSubsequence(UInt index, UInt number) const)
  AASequence seq1("DFPIANGER"), seq2("IAN"), seq3("DFPIANGER");
	TEST_EQUAL(seq2, seq1.getSubsequence(3, 3))
	TEST_EQUAL(seq3, seq1.getSubsequence(0, 9))
	TEST_EXCEPTION(Exception::IndexOverflow, seq1.getSubsequence(0, 10))
RESULT

CHECK(bool has(const Residue& residue) const)
  AASequence seq("DFPIANGER");
	TEST_EQUAL(seq.has(seq[0]), true)
	Residue res;
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

CHECK(bool hasSubsequence(const String& peptide) const)
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

CHECK(bool hasPrefix(const String& peptide) const)
  AASequence seq1("DFPIANGER");
	String seq2("DFP"), seq3("AIN"), seq4("BLUBB");
	TEST_EQUAL(seq1.hasPrefix(seq2), true)
	TEST_EQUAL(seq1.hasPrefix(seq3), false)
RESULT

CHECK(bool hasSuffix(const AASequence& peptide) const)
  AASequence seq1("DFPIANGER"), seq2("GER"), seq3("AIN");
  TEST_EQUAL(seq1.hasSuffix(seq2), true)
  TEST_EQUAL(seq1.hasSuffix(seq3), false) 
RESULT

CHECK(bool hasSuffix(const String& peptide) const)
  AASequence seq1("DFPIANGER");
  String seq2("GER"), seq3("AIN"), seq4("BLUBB");
  TEST_EQUAL(seq1.hasSuffix(seq2), true)
  TEST_EQUAL(seq1.hasSuffix(seq3), false)
RESULT

CHECK(ConstIterator begin() const)
	String result[] = { "D", "F", "P", "I", "A", "N", "G", "E", "R" };
	AASequence seq("DFPIANGER");
	UInt i = 0;
	for (AASequence::ConstIterator it = seq.begin(); it != seq.end(); ++it, ++i)
	{
		TEST_EQUAL((*it).getOneLetterCode(), result[i])
	}
RESULT

CHECK(ConstIterator end() const)
	NOT_TESTABLE
RESULT

CHECK(Iterator begin())
	String result[] = { "D", "F", "P", "I", "A", "N", "G", "E", "R" };
  AASequence seq("DFPIANGER");
  UInt i = 0;
  for (AASequence::ConstIterator it = seq.begin(); it != seq.end(); ++it, ++i)
  {
    TEST_EQUAL((*it).getOneLetterCode(), result[i])
  }
RESULT
		  
CHECK(Iterator end())
	NOT_TESTABLE
RESULT

//CHECK(friend std::ostream& operator << (std::ostream& os, const AASequence& peptide))
//  // TODO
//RESULT

//CHECK(friend std::istream& operator >> (std::istream& is, const AASequence& peptide))
//  // TODO
//RESULT

CHECK(String toString() const)
	AASequence seq1("DFPIANGER");
	AASequence seq2("(MOD:00051)DFPIANGER");
	AASequence seq3("DFPIAN(Deamidated)GER");

	TEST_EQUAL(seq1.isValid(), true)
	TEST_EQUAL(seq2.isValid(), true)
	TEST_EQUAL(seq3.isValid(), true)

	TEST_STRING_EQUAL(seq1.toString(), "DFPIANGER")
	TEST_STRING_EQUAL(seq2.toString(), "(MOD:00051)DFPIANGER")
	TEST_STRING_EQUAL(seq3.toString(), "DFPIAN(MOD:00565)GER")
RESULT

CHECK(String toUnmodifiedString() const)
	AASequence seq1("DFPIANGER");
  AASequence seq2("(MOD:00051)DFPIANGER");
  AASequence seq3("DFPIAN(Deamidated)GER");

  TEST_EQUAL(seq1.isValid(), true)
  TEST_EQUAL(seq2.isValid(), true)
  TEST_EQUAL(seq3.isValid(), true)

  TEST_STRING_EQUAL(seq1.toUnmodifiedString(), "DFPIANGER")
  TEST_STRING_EQUAL(seq2.toUnmodifiedString(), "DFPIANGER")
  TEST_STRING_EQUAL(seq3.toUnmodifiedString(), "DFPIANGER")
RESULT

CHECK(AASequence(const char *rhs))
	AASequence seq1('C');
	AASequence seq2('A');
	TEST_STRING_EQUAL(seq1.toString(), "C")
	TEST_STRING_EQUAL(seq2.toUnmodifiedString(), "A")
	AASequence seq3("CA");
	TEST_STRING_EQUAL((seq1 + seq2).toString(), seq3.toString())
RESULT

CHECK(void setModification(UInt index, const String &modification))
	AASequence seq1("ACDEFNK");
	seq1.setModification(5, "Deamidated");
	TEST_STRING_EQUAL(seq1[5].getModification(), "MOD:00565")
RESULT

CHECK(void setNTerminalModification(const String &modification))
	AASequence seq1("DFPIANGER");
	AASequence seq2("(MOD:00051)DFPIANGER");
	TEST_EQUAL(seq1 == seq2, false)
	seq1.setNTerminalModification("MOD:00051");
	TEST_EQUAL(seq1 == seq2, true)
	
	AASequence seq3("DABCDEF");
	AASequence seq4("(MOD:00051)DABCDEF");
	TEST_EQUAL(seq3 == seq4, false)
	TEST_EQUAL(seq3.isValid(), seq4.isValid())
	seq3.setNTerminalModification("MOD:00051");
	TEST_EQUAL(seq3.isModified(), true)
	TEST_EQUAL(seq4.isModified(), true)
	TEST_EQUAL(seq3 == seq4, true)
RESULT


CHECK(const String& getNTerminalModification() const)
	AASequence seq1("(MOD:00051)DFPIANGER");
	TEST_EQUAL(seq1.getNTerminalModification(), "MOD:00051")

	AASequence seq2("DFPIANGER");
	TEST_EQUAL(seq2.getNTerminalModification(), "")

RESULT


CHECK(void setCTerminalModification(const String &modification))
	AASequence seq1("DFPIANGER");
	AASequence seq2("DFPIANGER(ArgN)");

	TEST_EQUAL(seq1 == seq2, false)
	seq1.setCTerminalModification("ArgN");
	TEST_EQUAL(seq1 == seq2, true)

	AASequence seq3("DABCDER");
	AASequence seq4("DABCDER(ArgN)");
	TEST_EQUAL(seq3 == seq4, false)
	TEST_EQUAL(seq3.isValid(), seq4.isValid())
	seq3.setCTerminalModification("ArgN");
	TEST_EQUAL(seq3.isModified(), true)
	TEST_EQUAL(seq4.isModified(), true)
	TEST_EQUAL(seq3 == seq4, true)

	AASequence seq5("DABCDER(MOD:00177)");
	AASequence seq6("DABCDER(MOD:00177)(ArgN)");
	TEST_EQUAL(seq5.isModified(), false)
	TEST_EQUAL(seq6.isModified(), true)
	seq5.setCTerminalModification("ArgN");	
	TEST_EQUAL(seq5 == seq6, true)

	AASequence seq7("DFPIANGER(MOD:00177)");
	AASequence seq8("DFPIANGER(MOD:00177)(ArgN)");
	TEST_EQUAL(seq7.isModified(), true)
	TEST_EQUAL(seq8.isModified(), true)
	seq7.setCTerminalModification("ArgN");
	TEST_EQUAL(seq5 == seq6, true)
RESULT

CHECK(const String& getCTerminalModification() const)
  AASequence seq1("DFPIANGER(ArgN)");
  TEST_EQUAL(seq1.getCTerminalModification(), "MOD:00091")

  AASequence seq2("DFPIANGER");
  TEST_EQUAL(seq2.getCTerminalModification(), "")	
RESULT

CHECK(bool setStringSequence(const String &sequence))
	AASequence seq1("DFPIANGER");
	AASequence seq2("(MOD:00051)DFPIAK");

	AASequence seq3 = seq1;

	TEST_EQUAL(seq1 == seq3, true)
	
	seq3.setStringSequence("(MOD:00051)DFPIAK");
	TEST_EQUAL(seq2 == seq3, true)

	seq3.setStringSequence("DFPIANGER");
	TEST_EQUAL(seq1 == seq3, true)	
RESULT

CHECK(AASequence operator + (const char *rhs) const)
  AASequence seq1("DFPIANGER"), seq2("DFP");
	TEST_EQUAL(seq1, seq2 + "IANGER");
RESULT

CHECK(AASequence& operator += (const char *rhs))
	AASequence seq1("DFPIANGER"), seq2("DFP");
	seq2 += "IANGER";
	TEST_EQUAL(seq1, seq2)
RESULT


CHECK(bool isValid() const)
	AASequence seq1("(MOD:00051)DABCDEF");
	AASequence seq2("DABCDEF");
	AASequence seq3("(MOD:00051)DFPIANGER");
	AASequence seq4("DFPIANGER");

	TEST_EQUAL(seq1.isValid(), false)
	TEST_EQUAL(seq2.isValid(), false)
	TEST_EQUAL(seq3.isValid(), true)
	TEST_EQUAL(seq4.isValid(), true)
RESULT

CHECK(bool hasNTerminalModification() const)
	AASequence seq1("(MOD:00051)DABCDEF");
	AASequence seq2("DABCDEF");

	TEST_EQUAL(seq1.hasNTerminalModification(), true)
	TEST_EQUAL(seq2.hasNTerminalModification(), false)
	
	AASequence seq3("(MOD:00051)DFPIANGER");
	AASequence seq4("DFPIANGER");
	TEST_EQUAL(seq3.hasNTerminalModification(), true)
	TEST_EQUAL(seq4.hasNTerminalModification(), false)
RESULT
 
CHECK(bool hasCTerminalModification() const)
	AASequence seq1("DFPIANGER(ArgN)");
	AASequence seq2("DFPIANGER");
	TEST_EQUAL(seq1.hasCTerminalModification(), true)
	TEST_EQUAL(seq2.hasCTerminalModification(), false)
	seq1.setCTerminalModification("");
	TEST_EQUAL(seq1.hasCTerminalModification(), false)
RESULT

CHECK(bool isModified() const)
	AASequence seq1("DFPIANGER");
	TEST_EQUAL(seq1.isModified(), false);
	AASequence seq2(seq1);
	seq2.setNTerminalModification("MOD:09999");
	TEST_EQUAL(seq2.isModified(), true)

	AASequence seq3(seq1);
	seq3.setCTerminalModification("ArgN");
	TEST_EQUAL(seq3.isModified(), true);

	AASequence seq4("DFPIANGER(MOD:00177)");
	TEST_EQUAL(seq4.isModified(), true);
RESULT

CHECK(bool isModified(UInt index) const)
	AASequence seq4("DFPIAN(MOD:00565)GER");
	TEST_EQUAL(seq4.isModified(5), true)
	TEST_EQUAL(seq4.isModified(4), false)
RESULT

CHECK(bool operator<(const AASequence &rhs) const)
	AASequence seq1("DFPIANGER");
	AASequence seq2("DFBIANGER");
	TEST_EQUAL(seq2 < seq1, true)
	TEST_EQUAL(seq1 < seq2, false)
	AASequence seq3("DFPIANGFR");
	TEST_EQUAL(seq3 < seq1, false)
RESULT
 
CHECK(bool operator!=(const AASequence& rhs) const)
  AASequence seq1("(MOD:00051)DFPIANGER");
  AASequence seq2("DFPIANGER");
  TEST_EQUAL(seq2 != AASequence("DFPIANGER"), false)
  TEST_EQUAL(seq1 != AASequence("(MOD:00051)DFPIANGER"), false)

  AASequence seq3("DFPIANGER(MOD:00177)");
  AASequence seq4("DFPIANGER(ArgN)");
  TEST_EQUAL(seq3 != AASequence("DFPIANGER"), true)
  TEST_EQUAL(seq3 != AASequence("DFPIANGER(MOD:00177)"), false)
  TEST_EQUAL(seq4 != AASequence("DFPIANGER(ArgN)"), false)
  TEST_EQUAL(seq4 != AASequence("DFPIANGER"), true)

  AASequence seq5("DFBIANGER");
  TEST_EQUAL(seq5 != AASequence("DFPIANGER"), true)
  TEST_EQUAL(seq5 != AASequence("DFBIANGER"), false)
RESULT

CHECK(bool operator!=(const String& rhs) const)
  AASequence seq1("(MOD:00051)DFPIANGER");
  AASequence seq2("DFPIANGER");
  TEST_EQUAL(seq2 != String("DFPIANGER"), false)
  TEST_EQUAL(seq1 != String("(MOD:00051)DFPIANGER"), false)

  AASequence seq3("DFPIANGER(MOD:00177)");
  AASequence seq4("DFPIANGER(ArgN)");
  TEST_EQUAL(seq3 != String("DFPIANGER"), true)
  TEST_EQUAL(seq3 != String("DFPIANGER(MOD:00177)"), false)
  TEST_EQUAL(seq4 != String("DFPIANGER(ArgN)"), false)
  TEST_EQUAL(seq4 != String("DFPIANGER"), true)

	AASequence seq5("DFBIANGER");
	TEST_EQUAL(seq5 != String("DFPIANGER"), true)
	TEST_EQUAL(seq5 != String("DFBIANGER"), false)
RESULT

CHECK(bool operator!=(const char *rhs) const)
	AASequence seq1("(MOD:00051)DFPIANGER");
	AASequence seq2("DFPIANGER");
	TEST_EQUAL(seq2 != "DFPIANGER", false)
	TEST_EQUAL(seq1 != "(MOD:00051)DFPIANGER", false)

	AASequence seq3("DFPIANGER(MOD:00177)");
	AASequence seq4("DFPIANGER(ArgN)");
	TEST_EQUAL(seq3 != "DFPIANGER", true)
	TEST_EQUAL(seq3 != "DFPIANGER(MOD:00177)", false)
	TEST_EQUAL(seq4 != "DFPIANGER(ArgN)", false)
	TEST_EQUAL(seq4 != "DFPIANGER", true)

	AASequence seq5("DFBIANGER");
	TEST_EQUAL(seq5 != "DFPIANGER", true)
	TEST_EQUAL(seq5 != "DFBIANGER", false)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
