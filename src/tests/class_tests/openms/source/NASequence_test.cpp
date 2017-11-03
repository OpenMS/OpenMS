// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
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

NASequence* ptr = 0;
NASequence* null_ptr = 0;
RibonucleotideDB* db = RibonucleotideDB::getInstance();

START_SECTION((NASequence()=default))
{
  ptr = new NASequence();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->getFivePrimeMod(), (void*) nullptr);
  TEST_EQUAL(ptr->getThreePrimeMod(), (void*) nullptr);
  TEST_EQUAL(ptr->size(), 0);
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

START_SECTION((bool operator<(const NASequence& rhs) const ))
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
}
END_SECTION

START_SECTION((void setFivePrimeMod(const RibonucleotideChainEnd *r)))
{
  NASequence aaa = NASequence::fromString("AAA");
  TEST_EQUAL(aaa.hasFivePrimeMod(), false);
  aaa.setFivePrimeMod(db->getRibonucleotide("(pN)"));  // 5' phosphate
  TEST_EQUAL(aaa.hasFivePrimeMod(), true);
  TEST_EQUAL(aaa.getFivePrimeMod()->getCode(), "(pN)");
  TEST_STRING_EQUAL(aaa.toString(), "[(pN)]AAA");
}
END_SECTION

START_SECTION((const RibonucleotideChainEnd* getFivePrimeMod() const ))
{
  // tested via (const RibonucleotideChainEnd* getFivePrimeMod() const )
  TEST_EQUAL(true, true);
}
END_SECTION

START_SECTION((void setThreePrimeMod(const RibonucleotideChainEnd* r)))
{
  NASequence aaa = NASequence::fromString("AAA");
  TEST_EQUAL(aaa.hasThreePrimeMod(), false);
  aaa.setThreePrimeMod(db->getRibonucleotide("(pN)"));
  TEST_EQUAL(aaa.hasThreePrimeMod(), true);
  TEST_EQUAL(aaa.getThreePrimeMod()->getCode(), "(pN)");
  TEST_STRING_EQUAL(aaa.toString(), "AAA[(pN)]");
}
END_SECTION

START_SECTION((const RibonucleotideChainEnd* getThreePrimeMod() const ))
{
  // tested via (void setThreePrimeMod(const RibonucleotideChainEnd* r))
  NOT_TESTABLE
}
END_SECTION

START_SECTION((bool hasThreePrimeMod() const ))
{
  // tested via (void setThreePrimeMod(const RibonucleotideChainEnd* r))
  NOT_TESTABLE
}
END_SECTION

START_SECTION((double getMonoWeight(NASequence::NASFragmentType type = NASequence::Full, Int charge = 0) const))
{
  NASequence seq = NASequence::fromString("GGG");
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::AminusB, -1), 803.117);
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::WIon, -1), 1052.143);
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::YIon, -1), 972.177);
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::DIon, -1), 1034.133);
  TEST_REAL_SIMILAR(seq.getMonoWeight(NASequence::AminusB, -2), 802.117);
  NASequence seq_notSym = NASequence::fromString("GAU");
  TEST_REAL_SIMILAR(seq_notSym.getMonoWeight(NASequence::AminusB, -1), 787.122);


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

START_SECTION((EmpiricalFormula getFormula(NASequence::NASFragmentType type = NASequence::Full, Int charge = 0) const))
{
  NASequence seq = NASequence::fromString("GG");
  TEST_EQUAL(seq.getFormula(NASequence::Full, -1), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C10H12N5O5"));
  TEST_EQUAL(seq.getFormula(NASequence::Full, -2), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C10H11N5O5"));
  TEST_EQUAL(seq.getFormula(NASequence::WIon, -1), EmpiricalFormula("C20H25N10O15P2"));
  TEST_EQUAL(seq.getFormula(NASequence::XIon, -1), EmpiricalFormula("C20H25N10O14P2"));
  TEST_EQUAL(seq.getFormula(NASequence::YIon, -1), EmpiricalFormula("C10H12N5O6P") + EmpiricalFormula("C10H12N5O6"));
  TEST_EQUAL(seq.getFormula(NASequence::ZIon, -1), EmpiricalFormula("C20H24N10O11P"));
  TEST_EQUAL(seq.getFormula(NASequence::AIon, -1), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C10H10N5O4"));
  TEST_EQUAL(seq.getFormula(NASequence::BIon, -1), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C10H12N5O5"));
  TEST_EQUAL(seq.getFormula(NASequence::CIon, -1), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C10H13N5O7P"));
  TEST_EQUAL(seq.getFormula(NASequence::DIon, -1), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C10H11N5O7P"));
  TEST_EQUAL(seq.getFormula(NASequence::AminusB, -1), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C5H5O3"));

  seq = NASequence::fromString("GGG");
  TEST_EQUAL(seq.getFormula(NASequence::Full, -1), EmpiricalFormula("C30H36N15O19P2"));
  TEST_EQUAL(seq.getFormula(NASequence::FivePrime, -1), EmpiricalFormula("C30H36N15O19P2"));
  //TEST_EQUAL(seq.getFormula(NASequence::CTerminal, 1), EmpiricalFormula("C10H12N5O7P") + EmpiricalFormula("C10H12N5O7P"));
  //TEST_EQUAL(seq.getFormula(NASequence::Internal, 1), EmpiricalFormula("C10H10N5O6P") * 3 - EmpiricalFormula("H"));
  TEST_EQUAL(seq.getFormula(NASequence::WIon, -2), EmpiricalFormula("C20H25N10O15P2") + EmpiricalFormula("C10H11N5O7P"));
  TEST_EQUAL(seq.getFormula(NASequence::XIon, -2), EmpiricalFormula("C20H25N10O14P2") + EmpiricalFormula("C10H11N5O7P"));
  TEST_EQUAL(seq.getFormula(NASequence::YIon, -2), EmpiricalFormula("C10H12N5O6P") + EmpiricalFormula("C10H12N5O6")+EmpiricalFormula("C10H11N5O7P"));
  TEST_EQUAL(seq.getFormula(NASequence::ZIon, -2), EmpiricalFormula("C20H24N10O11P") + EmpiricalFormula("C10H11N5O7P"));
}
END_SECTION

START_SECTION((static NASequence fromString(const String &s)))
{
  NASequence seq = NASequence::fromString(String("CUA"));
  TEST_STRING_EQUAL(seq.toString(), "CUA");
}
END_SECTION

START_SECTION((static NASequence fromString(const char *s)))
{
  NASequence seq = NASequence::fromString("GG");
  TEST_STRING_EQUAL(seq.toString(), "GG");
}
END_SECTION

START_SECTION((string toString()))
{
  NASequence seq = NASequence::fromString("GG");
  TEST_STRING_EQUAL(seq.toString(), "GG");
}
END_SECTION


START_SECTION((NASequence getPrefix(Size length) const))
{
  // TODO
}
END_SECTION

START_SECTION((NASequence getSuffix(Size length) const))
{
  // TODO
}
END_SECTION

START_SECTION((Iterator begin()))
{
  // TODO
}
END_SECTION

START_SECTION((ConstIterator begin() const))
{
  // TODO
}
END_SECTION

START_SECTION((Iterator end()))
{
  // TODO
}
END_SECTION

START_SECTION((ConstIterator end() const))
{
  // TODO
}
END_SECTION

START_SECTION((ConstIterator cbegin() const))
{
  // TODO
}
END_SECTION

START_SECTION((ConstIterator cend() const))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] ConstIterator() = default))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] ConstIterator(const std::vector<const Ribonucleotide*>* vec_ptr, difference_type position)))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] ConstIterator(const ConstIterator& rhs)))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] ConstIterator(const NASequence::Iterator& rhs)))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] virtual ~ConstIterator()))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] const_reference operator*() const))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] const_pointer operator->() const))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] const ConstIterator operator+(difference_type diff) const))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] difference_type operator-(ConstIterator rhs) const))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] const ConstIterator operator-(difference_type diff) const))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] bool operator==(const ConstIterator &rhs) const))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] bool operator!=(const ConstIterator &rhs) const))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] ConstIterator& operator++()))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] ConstIterator& operator--()))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::ConstIterator] ConstIterator& operator=(const ConstIterator& rhs)))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] Iterator()=default))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] Iterator(std::vector<const Ribonucleotide*>* vec_ptr, difference_type position)))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] Iterator(const Iterator& rhs)))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] virtual ~Iterator()))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] const_reference operator*() const ))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] const_pointer operator->() const ))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] pointer operator->()))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] const Iterator operator+(difference_type diff) const ))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] difference_type operator-(Iterator rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] const Iterator operator-(difference_type diff) const ))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] bool operator==(const Iterator& rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] bool operator!=(const Iterator& rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] Iterator& operator++()))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] Iterator& operator--()))
{
  // TODO
}
END_SECTION

START_SECTION(([NASequence::Iterator] Iterator& operator=(const Iterator& rhs)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
