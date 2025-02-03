// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Timo Sachsenberg, Samuel Wein $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <iostream>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(NASequence, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

NASequence* ptr = nullptr;
NASequence* null_ptr = nullptr;
RibonucleotideDB* db = RibonucleotideDB::getInstance();

START_SECTION((NASequence()=default))
{
  ptr = new NASequence();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->getFivePrimeMod(), (void*) nullptr);
  TEST_EQUAL(ptr->getThreePrimeMod(), (void*) nullptr);
  TEST_EQUAL(ptr->size(), 0);
  delete ptr;
}
END_SECTION


START_SECTION((NASequence(const NASequence&) = default))
{
  // test Copy Constructor
  NASequence aaa = NASequence::fromString("AAA");
  NASequence aaa2(aaa);

  TEST_EQUAL(aaa.size(), 3);
  TEST_EQUAL(aaa2.size(), 3);
  TEST_EQUAL(aaa.operator==(aaa2), true);
}
END_SECTION

START_SECTION((NASequence(NASequence&&) = default))
{
  // test Move constructor
  NASequence aaa(NASequence::fromString("AAA"));
  TEST_EQUAL(aaa.size(), 3);
}
END_SECTION

START_SECTION((NASequence& operator=(const NASequence&)& = default))
{
  // test Copy Assignment
  NASequence aaa = NASequence::fromString("AAA");
  NASequence c = NASequence::fromString("C");
  c = aaa;
  TEST_EQUAL(aaa.size(), 3);
  TEST_EQUAL(c.size(), 3);
  TEST_EQUAL(aaa.operator==(c), true);
}
END_SECTION

START_SECTION((NASequence& operator=(NASequence&&)& = default))
{
  // test Move Assignment
  NASequence c = NASequence::fromString("C");
  c = NASequence::fromString("AAA");
  TEST_EQUAL(c.size(), 3);
}
END_SECTION

START_SECTION((NASequence(std::vector<const Ribonucleotide*> s, const RibonucleotideChainEnd* five_prime, const RibonucleotideChainEnd* three_prime)))
{
  NASequence aaa = NASequence::fromString("AAA");
  NASequence aaa2(aaa.getSequence(), nullptr, nullptr);
  TEST_EQUAL(aaa2.operator==(aaa), true);
}
END_SECTION

START_SECTION((virtual ~NASequence()=default))
{
  NASequence* n = new NASequence();
  delete(n);
}
END_SECTION

START_SECTION((bool operator==(const NASequence& rhs) const))
{
  NASequence aaa = NASequence::fromString("AAA");
  NASequence aaa2(aaa);
  TEST_EQUAL(aaa.operator==(aaa2), true);
}
END_SECTION

START_SECTION((bool operator<(const NASequence& rhs) const))
{
  NASequence aaa = NASequence::fromString("AAA");
  NASequence aaaa = NASequence::fromString("AAAA");
  NASequence cccc = NASequence::fromString("CCCC");
  TEST_EQUAL(aaa.operator<(aaaa), true);
  TEST_EQUAL(aaaa.operator<(aaa), false);
  TEST_EQUAL(aaaa.operator<(cccc), true);
}
END_SECTION

START_SECTION((void setSequence(const std::vector<const Ribonucleotide*>& s)))
{
  NASequence aaa = NASequence::fromString("AAA");
  NASequence cccc = NASequence::fromString("CCCC");
  aaa.setSequence(cccc.getSequence());
  TEST_EQUAL(aaa.operator==(cccc), true);
  TEST_EQUAL(aaa.size(), 4);
}
END_SECTION

START_SECTION((const std::vector<const Ribonucleotide*>& getSequence() const ))
{
  // tested via (void setSequence(const std::vector< const Ribonucleotide * > &s))
  TEST_EQUAL(true, true);
}
END_SECTION

START_SECTION((std::vector<const Ribonucleotide*>& getSequence()))
{
  // tested via (void setSequence(const std::vector< const Ribonucleotide * > &s))
  TEST_EQUAL(true, true);
}
END_SECTION

START_SECTION((void set(size_t index, const Ribonucleotide* r)))
{
  NASequence aaaa = NASequence::fromString("AAAA");
  NASequence cccc = NASequence::fromString("CCCC");
  aaaa.set(2, cccc.get(2));
  TEST_EQUAL(aaaa, NASequence::fromString("AACA"));
}
END_SECTION

START_SECTION((const Ribonucleotide* get(size_t index)))
{
  // tested via ((void set(size_t index, const Ribonucleotide *r)))
  TEST_EQUAL(true, true);
}
END_SECTION

START_SECTION((const Ribonucleotide*& operator[](size_t index)))
{
  NASequence aaa = NASequence::fromString("AAA");
  const NASequence ggg = NASequence::fromString("GGG");
  aaa[1] = ggg[2];
  TEST_EQUAL(aaa, NASequence::fromString("AGA"));
}
END_SECTION

START_SECTION((const Ribonucleotide*& operator[](size_t index) const ))
{
  NASequence aaa = NASequence::fromString("AAA");
  const NASequence ggg = NASequence::fromString("GGG");
  aaa[1] = ggg[2];
  TEST_EQUAL(aaa[1], ggg[0]);
}
END_SECTION

START_SECTION((bool empty() const ))
{
  NASequence aaa;
  TEST_EQUAL(aaa.empty(), true);
}
END_SECTION

START_SECTION((size_t size() const ))
{
  NASequence seq;
  TEST_EQUAL(seq.size(), 0);
  seq = NASequence::fromString("UGG");
  TEST_EQUAL(seq.size(), 3);
  // don't count terminal phosphate in sequence length:
  seq = NASequence::fromString("pUGG");
  TEST_EQUAL(seq.size(), 3);
  seq = NASequence::fromString("*UGG");
  TEST_EQUAL(seq.size(), 3);
  seq = NASequence::fromString("UGGp");
  TEST_EQUAL(seq.size(), 3);
  seq = NASequence::fromString("pUGGp");
  TEST_EQUAL(seq.size(), 3);
}
END_SECTION

START_SECTION((void clear()))
{
  NASequence aaa = NASequence::fromString("AAA");
  aaa.clear();
  TEST_EQUAL(aaa.empty(), true);
}
END_SECTION

START_SECTION((bool hasFivePrimeMod() const ))
{
  NASequence aaa = NASequence::fromString("AAA");
  TEST_EQUAL(aaa.hasFivePrimeMod(), false);
  aaa = NASequence::fromString("pAAA");
  TEST_EQUAL(aaa.hasFivePrimeMod(), true);
  aaa = NASequence::fromString("*AAA");
  TEST_EQUAL(aaa.hasFivePrimeMod(), true);
}
END_SECTION

START_SECTION((void setFivePrimeMod(const RibonucleotideChainEnd* r)))
{
  NASequence aaa = NASequence::fromString("AAA");
  TEST_EQUAL(aaa.hasFivePrimeMod(), false);
  aaa.setFivePrimeMod(db->getRibonucleotide("pN"));  // 5' phosphate
  TEST_EQUAL(aaa.hasFivePrimeMod(), true);
  TEST_EQUAL(aaa.getFivePrimeMod()->getCode(), "pN");
  TEST_STRING_EQUAL(aaa.toString(), "[pN]AAA");
}
END_SECTION

START_SECTION((const RibonucleotideChainEnd* getFivePrimeMod() const))
{
  // tested via (const RibonucleotideChainEnd* getFivePrimeMod() const )
  TEST_EQUAL(true, true);
}
END_SECTION

START_SECTION((void setThreePrimeMod(const RibonucleotideChainEnd* r)))
{
  NASequence aaa = NASequence::fromString("AAA");
  TEST_EQUAL(aaa.hasThreePrimeMod(), false);
  aaa.setThreePrimeMod(db->getRibonucleotide("N2'3'cp"));
  TEST_EQUAL(aaa.hasThreePrimeMod(), true);
  TEST_EQUAL(aaa.getThreePrimeMod()->getCode(), "N2'3'cp");
  TEST_STRING_EQUAL(aaa.toString(), "AAA[N2'3'cp]");
}
END_SECTION

START_SECTION((const RibonucleotideChainEnd* getThreePrimeMod() const))
{
  // tested via (void setThreePrimeMod(const RibonucleotideChainEnd* r))
  NOT_TESTABLE
}
END_SECTION

START_SECTION((bool hasThreePrimeMod() const))
{
  // tested via (void setThreePrimeMod(const RibonucleotideChainEnd* r))
  NOT_TESTABLE
}
END_SECTION


START_SECTION((EmpiricalFormula getFormula(NASequence::NASFragmentType type = NASequence::Full, Int charge = 0) const))
{
  NASequence seq = NASequence::fromString("GG");
  TEST_EQUAL(seq.getFormula(NASequence::Full, -1), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C10H12N5O5"));
  TEST_EQUAL(seq.getFormula(NASequence::Full, -2), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C10H11N5O5"));
  TEST_EQUAL(seq.getFormula(NASequence::WIon, -1), EmpiricalFormula("C20H25N10O15P2"));
  TEST_EQUAL(seq.getFormula(NASequence::XIon, -1), EmpiricalFormula("C20H23N10O14P2"));
  TEST_EQUAL(seq.getFormula(NASequence::YIon, -1), EmpiricalFormula("C10H12N5O6P") + EmpiricalFormula("C10H12N5O6"));
  TEST_EQUAL(seq.getFormula(NASequence::ZIon, -1), EmpiricalFormula("C20H22N10O11P"));
  TEST_EQUAL(seq.getFormula(NASequence::AIon, -1), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C10H10N5O4"));
  TEST_EQUAL(seq.getFormula(NASequence::BIon, -1), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C10H12N5O5"));
  TEST_EQUAL(seq.getFormula(NASequence::CIon, -1), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C10H11N5O7P"));
  TEST_EQUAL(seq.getFormula(NASequence::DIon, -1), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C10H13N5O8P"));
  TEST_EQUAL(seq.getFormula(NASequence::AminusB, -1), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C5H5O3"));
  // Thiol stuff
  NASequence thiolseq = NASequence::fromString("*[dA*][dA]");
  TEST_EQUAL(thiolseq.getFormula(NASequence::WIon, -1), EmpiricalFormula("C20H25N10O9P2S2"));
  TEST_EQUAL(thiolseq.getFormula(NASequence::YIon, -1), EmpiricalFormula("C20H24N10O7P1S1"));
  // Test for a-B ions
  thiolseq = NASequence::fromString("[dA][dA*]");
  TEST_EQUAL(thiolseq.getFormula(NASequence::AminusB, 0), EmpiricalFormula("C15H18N5O7P1"));
  thiolseq = NASequence::fromString("[dA][dA*][dA*]");
  TEST_EQUAL(thiolseq.getFormula(NASequence::AminusB, 0), EmpiricalFormula("C25H30N10O11P2S1"));
}
END_SECTION

START_SECTION((double getMonoWeight(NASequence::NASFragmentType type = NASequence::Full, Int charge = 0) const))
{
  // masses from Mongo-Oligo (http://mods.rna.albany.edu/masspec/Mongo-Oligo):
  NASequence seq = NASequence::fromString("GGG");
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::AminusB, -1), 803.117);
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::WIon, -1), 1052.143);
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::YIon, -1), 972.177);
  // Mongo-Oligo calls this ion "d-H20":
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::CIon, -1), 1034.133);
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::AminusB, -2), 802.117);
  NASequence seq_not_sym = NASequence::fromString("GAU");
  TEST_REAL_SIMILAR(seq_not_sym.getMonoWeight(NASequence::AminusB, -1), 787.122);

  seq = NASequence::fromString("AAUC");
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::AminusB, -1), 1077.1548);
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::CIon, -1), 1268.1644);
  seq = NASequence::fromString("AUCGp");
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::WIon, -1), 1382.1362);
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::YIon, -1), 1302.1698);

  seq = NASequence::fromString("[m1A]UCCACA");
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::AminusB, -1), 2006.2943);
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::CIon, -1), 2221.3151);
  seq = NASequence::fromString("UCCACAGp");
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::WIon, -1), 2321.2713);
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::YIon, -1), 2241.3049);

  // these masses were checked against external tools:
  seq = NASequence::fromString("pAAUCCAUGp");
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::Full, 0), 2652.312);
  seq = NASequence::fromString("ACCAAAGp");
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::Full, 0), 2289.348);
  seq = NASequence::fromString("AUUCACCC");
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::Full, 0), 2428.362);

  // with charge (negative!):
  seq = NASequence::fromString("AAU[m5C]Gp");
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::Full, -2), 1644.228);
}
END_SECTION

START_SECTION((double getAverageWeight(NASequence::NASFragmentType type = NASequence::Full, Int charge = 0) const))
{
  // data from RNAModMapper publication (Yu et al., Anal. Chem. 2017), Fig. 4:
  NASequence seq = NASequence::fromString("A[ms2i6A]AACCGp");
  TEST_REAL_SIMILAR(seq.getAverageWeight(NASequence::Full, -2) / 2, 1201.3);
  seq = NASequence::fromString("A[ms2i6A]AACC");
  TEST_REAL_SIMILAR(seq.getAverageWeight(NASequence::AminusB, -1), 1848.587 + 0.734);
  TEST_REAL_SIMILAR(seq.getAverageWeight(NASequence::CIon, -2) / 2, 1020.023 - 0.324);
  seq = NASequence::fromString("[ms2i6A]AACCGp");
  TEST_REAL_SIMILAR(seq.getAverageWeight(NASequence::WIon, -2) / 2, 1076.045 + 0.651);
  TEST_REAL_SIMILAR(seq.getAverageWeight(NASequence::YIon, -2) / 2, 1036.459 + 0.247);
}
END_SECTION

START_SECTION((static NASequence fromString(const String& s)))
{
  NASequence seq = NASequence::fromString(String("CUA"));
  TEST_STRING_EQUAL(seq.toString(), "CUA");
}
END_SECTION

START_SECTION((static NASequence fromString(const char* s)))
{
  NASequence seq = NASequence::fromString("GG");
  TEST_STRING_EQUAL(seq.toString(), "GG");
}
END_SECTION

START_SECTION((string toString()))
{
  NASequence seq = NASequence::fromString("GG");
  TEST_STRING_EQUAL(seq.toString(), "GG");
    // Test that all the entries in the Mod file work can be parsed from strings
  for(auto it = db->begin(); it != db->end(); ++it)
  {
    if ((*it)->getCode() != "c") // handle the phosphates and cyclophosphates separately, we pass on them for now
    {
      continue;
    } 
    else if((*it)->getCode() != "p")
    {
      continue;
    }
    else
    {
      //check that we can get the code from the entry, then convert it to a string and back
      auto ribocode = NASequence().fromString("[" + (*it)->getCode() + "]"); // get the code, which may or may not have brackets
      String asString = ribocode.toString();
      if (asString.hasPrefix("[")) // if we were multicharacter toString adds brackets
      {
        asString = asString.substr(1, asString.size() - 2);
      }
      TEST_STRING_EQUAL(asString, (*it)->getCode());
    }
  }
}
END_SECTION


START_SECTION((NASequence getPrefix(Size length) const))
{
  NASequence seq = NASequence::fromString("A[ms2i6A]AACCGp");
  NASequence seq2 = NASequence::fromString("A[ms2i6A]");
  NASequence seq3 = NASequence::fromString("AAACCG");
  NASequence seq4 = NASequence::fromString("AAA");
  TEST_EQUAL(seq.getPrefix(2),seq2);
  TEST_EQUAL(seq3.getPrefix(3),seq4);
  TEST_NOT_EQUAL(seq.getPrefix(3),seq2);
  TEST_NOT_EQUAL(seq.getPrefix(3),seq4);
  TEST_EXCEPTION(Exception::IndexOverflow, seq.getPrefix(10));
}
END_SECTION

START_SECTION((NASequence getSuffix(Size length) const))
{
  NASequence seq = NASequence::fromString("A[ms2i6A]AACAGp");
  NASequence seq2 = NASequence::fromString("[ms2i6A]AACAGp");
  NASequence seq3 = NASequence::fromString("AAACAG");
  NASequence seq4 = NASequence::fromString("CAG");
  NASequence seq5 = NASequence::fromString("[C*]AG");
  NASequence seq6 = NASequence::fromString("*AG");
  TEST_EQUAL(seq.getSuffix(6),seq2);
  TEST_EQUAL(seq3.getSuffix(3),seq4);
  TEST_NOT_EQUAL(seq.getSuffix(3),seq2);
  TEST_NOT_EQUAL(seq.getSuffix(3),seq4);
  TEST_EXCEPTION(Exception::IndexOverflow, seq.getSuffix(10));
  // Check that we handle the weird case of thiol
  TEST_STRING_EQUAL(seq5.getSuffix(2).toString(), seq6.toString() );
}
END_SECTION

START_SECTION((NASequence getSubsequence(Size start, Size length) const))
{
  NASequence seq = NASequence::fromString("pAUCGp");
  NASequence seq5 = NASequence::fromString("[C*]AG");
  TEST_STRING_EQUAL(seq.getSubsequence().toString(), "pAUCGp");
  TEST_STRING_EQUAL(seq.getSubsequence(1).toString(), "UCGp");
  TEST_STRING_EQUAL(seq.getSubsequence(0, 2).toString(), "pAU");
  TEST_STRING_EQUAL(seq.getSubsequence(2, 1).toString(), "C");
  //thiol stuff
  TEST_STRING_EQUAL(seq5.getSubsequence(1).toString(), "*AG");
}
END_SECTION

START_SECTION((Iterator begin()))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=0;
  for (NASequence::Iterator it = seq.begin(); it != seq.end(); ++it, ++i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION((ConstIterator begin() const))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=0;
  for (NASequence::ConstIterator it = seq.begin(); it != seq.end(); ++it, ++i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION((Iterator end()))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=0;
  for (NASequence::Iterator it = seq.begin(); it != seq.end(); ++it, ++i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION((ConstIterator end() const))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=0;
  for (NASequence::ConstIterator it = seq.begin(); it != seq.end(); ++it, ++i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION((ConstIterator cbegin() const))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=0;
  for (NASequence::ConstIterator it = seq.cbegin(); it != seq.cend(); ++it, ++i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION((ConstIterator cend() const))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=0;
  for (NASequence::ConstIterator it = seq.cbegin(); it != seq.cend(); ++it, ++i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] ConstIterator() = default))
{
  NASequence::ConstIterator iter = NASequence::ConstIterator(); // fails if it segfault
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] ConstIterator(const std::vector<const Ribonucleotide*>* vec_ptr, difference_type position)))
{
  NASequence seq = NASequence::fromString("AUCG");
  NASequence::ConstIterator it = seq.cbegin();
  TEST_EQUAL((*(it+2)).getCode(), "C");
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] ConstIterator(const ConstIterator& rhs)))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] ConstIterator(const NASequence::Iterator& rhs)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] virtual ~ConstIterator()))
{
NASequence::ConstIterator* ptr = new NASequence::ConstIterator();
delete ptr;
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] const_reference operator*() const))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=0;
  for (NASequence::ConstIterator it = seq.cbegin(); it != seq.cend(); ++it, ++i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] const_pointer operator->() const))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] const ConstIterator operator+(difference_type diff) const))
{
  NASequence seq = NASequence::fromString("AUCG");
  NASequence::ConstIterator it = seq.begin();
  TEST_EQUAL((*(it+2)).getCode(), "C");
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] difference_type operator-(ConstIterator rhs) const))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] const ConstIterator operator-(difference_type diff) const))
{
  NASequence seq = NASequence::fromString("AUCG");
  NASequence::ConstIterator it = seq.end();
  TEST_EQUAL((*(it-2)).getCode(), "C");
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] bool operator==(const ConstIterator &rhs) const))
{
  NASequence seq = NASequence::fromString("AUCG");
  NASequence::ConstIterator it = seq.end();
  TEST_EQUAL((it-4 == seq.begin()), true);
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] bool operator!=(const ConstIterator &rhs) const))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=0;
  for (NASequence::ConstIterator it = seq.cbegin(); it != seq.cend(); ++it, ++i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] ConstIterator& operator++()))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=0;
  for (NASequence::ConstIterator it = seq.cbegin(); it != seq.cend(); ++it, ++i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] ConstIterator& operator--()))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=3;
  for (NASequence::ConstIterator it = seq.end()-1; it != seq.begin(); --it, --i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] ConstIterator& operator=(const ConstIterator& rhs)))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=0;
  for (NASequence::ConstIterator it = seq.cbegin(); it != seq.cend(); ++it, ++i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION(([NASequence::Iterator] Iterator()=default))
{
  NASequence::Iterator iter= NASequence::Iterator();
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([NASequence::Iterator] Iterator(std::vector<const Ribonucleotide*>* vec_ptr, difference_type position)))
{
  NASequence seq = NASequence::fromString("AUCG");
  NASequence::Iterator it = seq.begin();
  TEST_EQUAL((*(it+2)).getCode(), "C");
}
END_SECTION

START_SECTION(([NASequence::Iterator] Iterator(const Iterator& rhs)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([NASequence::Iterator] virtual ~Iterator()))
{
  NASequence::Iterator* ptr = new NASequence::Iterator();
  delete ptr;
}
END_SECTION

START_SECTION(([NASequence::Iterator] const_reference operator*() const))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=0;
  for (NASequence::Iterator it = seq.begin(); it != seq.end(); ++it, ++i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION(([NASequence::Iterator] const_pointer operator->() const))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([NASequence::Iterator] pointer operator->()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([NASequence::Iterator] const Iterator operator+(difference_type diff) const))
{
  NASequence seq = NASequence::fromString("AUCG");
  NASequence::Iterator it = seq.begin();
  TEST_EQUAL((*(it+2)).getCode(), "C");
}
END_SECTION

START_SECTION(([NASequence::Iterator] difference_type operator-(Iterator rhs) const))
{
  //TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] const Iterator operator-(difference_type diff) const))
{
  NASequence seq = NASequence::fromString("AUCG");
  NASequence::Iterator it = seq.end();
  TEST_EQUAL((*(it-2)).getCode(), "C");
}
END_SECTION

START_SECTION(([NASequence::Iterator] bool operator==(const Iterator& rhs) const))
{
  NASequence seq = NASequence::fromString("AUCG");
  NASequence::Iterator it = seq.end();
  TEST_EQUAL((it-4 == seq.begin()), true);
}
END_SECTION

START_SECTION(([NASequence::Iterator] bool operator!=(const Iterator& rhs) const))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=0;
  for (NASequence::Iterator it = seq.begin(); it != seq.end(); ++it, ++i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION(([NASequence::Iterator] Iterator& operator++()))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=0;
  for (NASequence::Iterator it = seq.begin(); it != seq.end(); ++it, ++i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION(([NASequence::Iterator] Iterator& operator--()))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=3;
  for (NASequence::Iterator it = seq.end()-1; it != seq.begin(); --it, --i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION

START_SECTION(([NASequence::Iterator] Iterator& operator=(const Iterator& rhs)))
{
  String result[] = {"A","U","C","G"};
  NASequence seq = NASequence::fromString("AUCG");
  Size i=0;
  for (NASequence::Iterator it = seq.begin(); it != seq.end(); ++it, ++i)
  {
    TEST_EQUAL((*it).getCode(), result[i]);
  }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
