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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Tobias Rausch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/ListUtils.h>


using namespace OpenMS;
using namespace std;

///////////////////////////
///////////////////////////

void setDB(const StringList& in, AhoCorasickAmbiguous::PeptideDB& out)
{
  clear(out);
  for (int i = 0; i< in.size(); ++i) 
  {
    seqan::appendValue(out, in[i].c_str());
  }
}

void compareHits(int line, String protein, String expected_s, StringList observed, AhoCorasickAmbiguous::PeptideDB& pep_db)
{
  std::cout << "results of test line " << line << " for protein " << protein << ":\n";
  StringList expected = ListUtils::create<String>(expected_s.removeWhitespaces(), ',');
  std::sort(expected.begin(), expected.end());
  std::sort(observed.begin(), observed.end());
  TEST_EQUAL(observed.size(), expected.size()) // results should have same number of entries
  if (expected.size() == observed.size())
  {
    for (int i = 0; i < expected.size(); ++i)
    {
      expected[i] = expected[i].toUpper();
      std::cout << "hit " << i << ": " << expected[i] << " <> " << observed[i] << "\n";
      TEST_EQUAL(observed[i], expected[i])
    }
  }
  else
  {
    std::cout << "Results differ in number of hits:\n  expected:\n    " << ListUtils::concatenate(expected, "\n    ") << "  \nobserved:\n    " << ListUtils::concatenate(observed, "\n    ") << "\n";
  }
}

START_TEST(AhoCorasickAmbiguous, "$Id$")

AhoCorasickAmbiguous* ptr = 0;
AhoCorasickAmbiguous* nullPointer = 0;
START_SECTION(AhoCorasickAmbiguous())
{
  ptr = new AhoCorasickAmbiguous("test");
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~AhoCorasickAmbiguous())
{
  delete ptr;
}
END_SECTION

AhoCorasickAmbiguous::FuzzyACPattern pattern;
AhoCorasickAmbiguous::PeptideDB pep_db;

START_SECTION(static void initPattern(const PeptideDB& pep_db, const int aaa_max, FuzzyACPattern& pattern))
  setDB(ListUtils::create<String>("withB", ','), pep_db); // ambiguous char in peptide DB not allowed
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue, AhoCorasickAmbiguous::initPattern(pep_db, 2, pattern), "The value 'WITHB' was used but is not valid! Input peptide to FuzzyAC must NOT contain ambiguous amino acids (B/J/Z/X)!")
  setDB(ListUtils::create<String>("withJ", ','), pep_db); // unknown chars are converted to 'X'; ambiguous char in peptide DB not allowed
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue, AhoCorasickAmbiguous::initPattern(pep_db, 2, pattern), "The value 'WITHJ' was used but is not valid! Input peptide to FuzzyAC must NOT contain ambiguous amino acids (B/J/Z/X)!")
  setDB(ListUtils::create<String>("withZ", ','), pep_db); // unknown chars are converted to 'X'; ambiguous char in peptide DB not allowed
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue, AhoCorasickAmbiguous::initPattern(pep_db, 2, pattern), "The value 'WITHZ' was used but is not valid! Input peptide to FuzzyAC must NOT contain ambiguous amino acids (B/J/Z/X)!")
  setDB(ListUtils::create<String>("withX", ','), pep_db); // unknown chars are converted to 'X'; ambiguous char in peptide DB not allowed
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue, AhoCorasickAmbiguous::initPattern(pep_db, 2, pattern), "The value 'WITHX' was used but is not valid! Input peptide to FuzzyAC must NOT contain ambiguous amino acids (B/J/Z/X)!")

  setDB(ListUtils::create<String>("acd,adc,cad,cda,dac,dca", ','), pep_db);
  AhoCorasickAmbiguous::initPattern(pep_db, 2, pattern);
  TEST_EQUAL((UInt)pattern.max_ambAA, 2)
  AhoCorasickAmbiguous::initPattern(pep_db, 3, pattern);
  TEST_EQUAL((UInt)pattern.max_ambAA, 3)

  TEST_EQUAL(seqan::numVertices(pattern.data_graph), 16); // 1 root,5x3 subtrees

END_SECTION

START_SECTION(AhoCorasickAmbiguous(const String& protein_sequence))
  AhoCorasickAmbiguous fuzzyAC("XXX");

  AhoCorasickAmbiguous fuzzyAC2("BXZU");

  AhoCorasickAmbiguous fuzzyAC3("BXZU");

END_SECTION

START_SECTION(bool findNext(const FuzzyACPattern& pattern))

  //
  // Note: this only finds infixes. we do not care about trypticity at this level!
  //
  String prot;
  String expected;
  AhoCorasickAmbiguous fuzzyAC;
  std::vector<String> observed;
  /////////////////////////
  // "acd,adc,cad,cda,dac,dca"
  /////////////////////////
  AhoCorasickAmbiguous::initPattern(pep_db, 0, pattern);
  observed.clear();
  prot = "acdIadcIcadIcdaIdacIdca";
  fuzzyAC.setProtein(prot);
  expected = "acd@0,  adc@4,  cad@8,  cda@12,  dac@16,  dca@20"; // all six hits, found without spawning(ambAA)
  while (fuzzyAC.findNext(pattern))
  {
    observed.push_back(String(pep_db[fuzzyAC.getHitDBIndex()].data_begin, pep_db[fuzzyAC.getHitDBIndex()].data_end) + "@" + fuzzyAC.getHitProteinPosition());
  }
  compareHits(__LINE__, "XXX", expected, observed, pep_db);
  ///
  /// same, but with ambAA's allowed (but not used)
  ///
  AhoCorasickAmbiguous::initPattern(pep_db, 3, pattern);
  observed.clear();
  prot = "acdIadcIcadIcdaIdacIdca";
  fuzzyAC.setProtein(prot);
  expected = "acd@0,  adc@4,  cad@8,  cda@12,  dac@16,  dca@20"; // all six hits, found without spawning(ambAA)
  while (fuzzyAC.findNext(pattern))
  {
    observed.push_back(String(pep_db[fuzzyAC.getHitDBIndex()].data_begin, pep_db[fuzzyAC.getHitDBIndex()].data_end) + "@" + fuzzyAC.getHitProteinPosition());
  }
  compareHits(__LINE__, "XXX", expected, observed, pep_db);
  ///
  /// all ambAA's
  ///
  observed.clear();
  prot = "XXX";
  fuzzyAC.setProtein(prot);
  expected = "dac@0,  cad@0,  cda@0,  dca@0,  adc@0,  acd@0"; // all six hits, found at first position
  while (fuzzyAC.findNext(pattern))
  {
    observed.push_back(String(pep_db[fuzzyAC.getHitDBIndex()].data_begin, pep_db[fuzzyAC.getHitDBIndex()].data_end) + "@" + fuzzyAC.getHitProteinPosition());
  }
  compareHits(__LINE__, "XXX", expected, observed, pep_db);
  ///
  /// with prefix
  ///
  observed.clear();
  prot = "aXXX";
  fuzzyAC.setProtein(prot);
  expected = "acd@0,  adc@0,  dac@1,  cad@1,  cda@1,  dca@1,  adc@1,  acd@1"; // 2 hits of aXX at first pos; all six hits, found at second position
  while (fuzzyAC.findNext(pattern))
  {
    observed.push_back(String(pep_db[fuzzyAC.getHitDBIndex()].data_begin, pep_db[fuzzyAC.getHitDBIndex()].data_end) + "@" + fuzzyAC.getHitProteinPosition());
  }
  compareHits(__LINE__, prot, expected, observed, pep_db);
  ///
  /// with preifx and B instead of X
  ///
  observed.clear();
  prot = "aXBX"; // B = D|N,  Z = E|Q
  fuzzyAC.setProtein(prot);
  expected = "acd@0,     cda@1, adc@1"; //  
  while (fuzzyAC.findNext(pattern))
  {
    observed.push_back(String(pep_db[fuzzyAC.getHitDBIndex()].data_begin, pep_db[fuzzyAC.getHitDBIndex()].data_end) + "@" + fuzzyAC.getHitProteinPosition());
  }
  compareHits(__LINE__, prot, expected, observed, pep_db);
  ///
  /// test with two ambAA's: nothing should be found
  ///
  AhoCorasickAmbiguous::initPattern(pep_db, 2, pattern);
  fuzzyAC.setProtein("XXX");
  TEST_EQUAL(fuzzyAC.findNext(pattern), false);
  ///
  /// only two hits (due to ambAA==2)
  ///
  observed.clear();
  prot = "aXXX";
  fuzzyAC.setProtein(prot);
  expected = "acd@0,  adc@0    "; // 2 hits of aXX at first pos; nothing at second pos (since that requires three AAA)
  while (fuzzyAC.findNext(pattern))
  {
    observed.push_back(String(pep_db[fuzzyAC.getHitDBIndex()].data_begin, pep_db[fuzzyAC.getHitDBIndex()].data_end) + "@" + fuzzyAC.getHitProteinPosition());
  }
  compareHits(__LINE__, prot, expected, observed, pep_db);
  ///
  /// with suffix
  ///
  observed.clear();
  prot = "XXXc";
  fuzzyAC.setProtein(prot);
  expected = "adc@1,  dac@1    "; // 2 hits of XXc at second pos; nothing at first pos (since that requires three AAA)
  while (fuzzyAC.findNext(pattern))
  {
    observed.push_back(String(pep_db[fuzzyAC.getHitDBIndex()].data_begin, pep_db[fuzzyAC.getHitDBIndex()].data_end) + "@" + fuzzyAC.getHitProteinPosition());
  }
  compareHits(__LINE__, prot, expected, observed, pep_db);

  ///
  ///  new peptide DB
  ///
  setDB(ListUtils::create<String>("eq,nd,llll", ','), pep_db);
  AhoCorasickAmbiguous::initPattern(pep_db, 2, pattern);
  ///
  /// hits across the protein
  ///
  observed.clear();
  prot = "aXXaBBkkZZlllllk";  // B = D|N,  Z = E|Q
  fuzzyAC.setProtein(prot);
  expected = "nd@1, nd@4, eq@1, eq@8, llll@10, llll@11   "; // both match XX@1, eq matches ZZ, nd matches BB
  while (fuzzyAC.findNext(pattern))
  {
    observed.push_back(String(pep_db[fuzzyAC.getHitDBIndex()].data_begin, pep_db[fuzzyAC.getHitDBIndex()].data_end) + "@" + fuzzyAC.getHitProteinPosition());
  }
  compareHits(__LINE__, prot, expected, observed, pep_db);


END_SECTION

START_SECTION(Size getHitDBIndex())
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(Size getHitProteinPosition())
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION([EXTRA]template<typename T> inline void _getSpawnRange(const AAcid c, T& idxFirst, T& idxLast))
  
  // test that our AAcid translation table is correct
  for (char c = 'A'; c <= 'Z'; ++c)
  {
    TEST_EQUAL(char(seqan::AAcid(c)), c)
  }
  for (char c = 'a'; c <= 'z'; ++c)
  {
    TEST_EQUAL(char(seqan::AAcid(c)), String(c).toUpper()[0])
  }

  // check correct range
  Int32 idx_first, idx_last;
  seqan::_getSpawnRange(seqan::AAcid('B'), idx_first, idx_last);
  TEST_EQUAL(seqan::AAcid(idx_first), seqan::AAcid('D'))
  TEST_EQUAL(seqan::AAcid(idx_last), seqan::AAcid('N'))
  TEST_EQUAL(idx_last - idx_first, 1)

  seqan::_getSpawnRange(seqan::AAcid('J'), idx_first, idx_last);
  TEST_EQUAL(seqan::AAcid(idx_first), seqan::AAcid('I'))
  TEST_EQUAL(seqan::AAcid(idx_last), seqan::AAcid('L'))
  TEST_EQUAL(idx_last - idx_first, 1)

  seqan::_getSpawnRange(seqan::AAcid('Z'), idx_first, idx_last);
  TEST_EQUAL(seqan::AAcid(idx_first), seqan::AAcid('E'))
  TEST_EQUAL(seqan::AAcid(idx_last), seqan::AAcid('Q'))
  TEST_EQUAL(idx_last - idx_first, 1)
  
  seqan::_getSpawnRange(seqan::AAcid('X'), idx_first, idx_last);
  TEST_EQUAL(seqan::AAcid(idx_first), seqan::AAcid('A'))
  TEST_EQUAL(seqan::AAcid(idx_last), seqan::AAcid('V'))
  TEST_EQUAL(idx_last - idx_first, 21) // 21 unambiguous AA's
  
  // neighbours
  TEST_EQUAL(ordValue(seqan::AAcid('B')) - ordValue(seqan::AAcid('J')), -1)
  TEST_EQUAL(ordValue(seqan::AAcid('J')) - ordValue(seqan::AAcid('Z')), -1)
  TEST_EQUAL(ordValue(seqan::AAcid('Z')) - ordValue(seqan::AAcid('X')), -1)
END_SECTION

START_SECTION([EXTRA]inline bool isAmbiguous(AAcid c))
{
  String amb = "BJZX"; // all amb AA's
  for (char c = 'A'; c <= 'Z'; ++c) // test all characters from A..Z
  {
    TEST_EQUAL(seqan::isAmbiguous(seqan::AAcid(c)), amb.has(c));
  }
  amb = "bjzx"; // all amb AA's
  for (char c = 'a'; c <= 'z'; ++c) // test all characters from A..Z
  {
    TEST_EQUAL(seqan::isAmbiguous(seqan::AAcid(c)), amb.has(c));
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
