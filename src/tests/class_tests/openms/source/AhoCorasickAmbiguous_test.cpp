// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>
///////////////////////////

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <array>
#include <cassert>
#include <string_view>

using namespace OpenMS;
using namespace std;

///////////////////////////
///////////////////////////


void compareHits(int line, const string& protein, String expected_s, StringList& observed)
{
  std::cout << "results of test line " << line << " for protein " << protein << ":\n";
  //expected_s.toUpper();
  StringList expected = ListUtils::create<String>(expected_s.removeWhitespaces(), ',');
  std::sort(expected.begin(), expected.end());
  std::sort(observed.begin(), observed.end());
  TEST_EQUAL(observed.size(), expected.size()) // results should have same number of entries
  if (expected.size() == observed.size())
  {
    for (size_t i = 0; i < expected.size(); ++i)
    {
      //expected[i] = expected[i].toUpper();
      std::cout << "hit " << i << ": " << expected[i] << " <> " << observed[i] << "\n";
      TEST_EQUAL(expected[i], observed[i])
      if (expected[i] != observed[i])
      {
        std::cout << "difference!" << expected[i] << observed[i] << '\n';
      }
    }
  }
  else
  {
    std::cout << "Results differ in number of hits:\n  expected:\n    " << ListUtils::concatenate(expected, "\n    ") << "  \nobserved:\n    " << ListUtils::concatenate(observed, "\n    ") << "\n";
  }
}

void testCase(const ACTrie& t, const string& protein, const string& expected, vector<string>& needles, int line)
{
  std::vector<String> observed;
  ACTrieState state;
  state.setQuery(protein);
  while (t.nextHits(state))
  {
    for (auto& hit : state.hits)
    {
      observed.push_back(needles[hit.needle_index] + "@" + String(hit.query_pos));
    }
  }
  compareHits(line, protein, expected, observed);
}


template<int SIZE>
void checkAAIterator(const std::array<AA, SIZE>& aa_array, const std::array<size_t, SIZE>& pos_array, ACTrieState& state)
{
  size_t i = 0;
  for (AA aa = state.nextValidAA(); aa.isValid(); aa = state.nextValidAA(), ++i)
  {
    TEST_EQUAL(aa == aa_array[i], true);
    TEST_EQUAL(state.textPos(), pos_array[i]);
  }
  TEST_EQUAL(aa_array.size(), i)
}


START_TEST(AhoCorasickAmbiguous, "$Id$")

ACTrie* ptr = 0;
ACTrie* nullPointer = 0;
START_SECTION(ACTrie(uint32_t max_aaa = 0, uint32_t max_mm = 0))
{
  ptr = new ACTrie(1, 0);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ACTrie())
{
  delete ptr;
}
END_SECTION

START_SECTION(void addNeedle(const std::string& needle))
  ACTrie t(1,2);
  t.addNeedle("WITHV"); // normal AA's are allowed
  t.addNeedle("WITHB"); // ambiguous char 'B' is allowed
  t.addNeedle("WITHJ"); // ambiguous char 'J' is allowed
  t.addNeedle("WITHZ"); // ambiguous char 'Z' is allowed
  t.addNeedle("WITHX"); // ambiguous char 'X' is allowed
  TEST_EXCEPTION(Exception::InvalidValue, t.addNeedle("WITH$")) 
  TEST_EXCEPTION(Exception::InvalidValue, t.addNeedle("WITH*")) 

  TEST_EQUAL(t.getNeedleCount(), 5)
  TEST_EQUAL(t.getMaxAAACount(), 1)
  TEST_EQUAL(t.getMaxMMCount(), 2)
END_SECTION

START_SECTION(void addNeedles(const std::vector<std::string>& needle))
  ACTrie t(1, 2);
  t.addNeedle("WITHV"); // normal AA's are allowed
  t.addNeedle("WITHB"); // ambiguous char 'B' is allowed
  t.addNeedle("WITHJ"); // ambiguous char 'J' is allowed
  TEST_EXCEPTION(Exception::InvalidValue, t.addNeedles({"WITHZ", "WITHX","WITH$"}))

  TEST_EQUAL(t.getNeedleCount(), 5)
END_SECTION

START_SECTION(void compressTrie())
  NOT_TESTABLE // needs context...            
END_SECTION


START_SECTION(size_t getNeedleCount() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setMaxAAACount(const uint32_t max_aaa))
  NOT_TESTABLE // tested below
END_SECTION

START_SECTION(uint32_t getMaxAAACount() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setMaxMMCount(const uint32_t max_mm))
  NOT_TESTABLE // tested below
END_SECTION

START_SECTION(uint32_t getMaxMMCount() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(bool nextHits(ACTrieState& state) const)
{
  //
  // Note: we do not care about trypticity at this level!
  //
  ACTrie t(0, 0);
  std::vector<string> needles = {"acd", "adc", "cad", "cda", "dac", "dca"};
  t.addNeedlesAndCompress(needles);

  /////////////////////////
  // "acd,adc,cad,cda,dac,dca"
  /////////////////////////

  // all six hits, found without spawning(ambAA))
  testCase(t, "acdIadcIcadIcdaIdacIdca", "acd@0,  adc@4,  cad@8,  cda@12,  dac@16,  dca@20", needles, __LINE__);

  ///
  /// same, but with ambAA's allowed (but not used)
  ///
  t.setMaxAAACount(3);
  // all six hits, found without spawning(ambAA)
  testCase(t, "acdIadcIcadIcdaIdacIdca", "acd@0,  adc@4,  cad@8,  cda@12,  dac@16,  dca@20", needles, __LINE__);
  
  ///
  /// all ambAA's
  ///
  // all six hits, found at first position
  testCase(t, "XXX", "dac@0,  cad@0,  cda@0,  dca@0,  adc@0,  acd@0", needles, __LINE__);
 
  ///
  /// with prefix
  ///
  // 2 hits of aXX at first pos; all six hits, found at second position
  testCase(t, "aXXX", "acd@0, adc@0, dac@1, cad@1, cda@1, dca@1, adc@1, acd@1", needles, __LINE__);
  
  ///
  /// with prefix and B instead of X
  ///
  // B = D|N;  Z = E|Q
  testCase(t, "aXBX", "acd@0, cda@1, adc@1", needles, __LINE__);

  ///
  /// test with two ambAA's: nothing should be found
  ///
  t.setMaxAAACount(2);
  testCase(t, "XXX", "", needles, __LINE__);

  ///
  /// only two hits (due to ambAA==2)
  ///
  t.setMaxAAACount(2);
  // 2 hits of aXX at first pos; nothing at second pos (since that requires three AAA)
  testCase(t, "aXXX", "acd@0,  adc@0", needles, __LINE__);

  ///
  /// with suffix
  ///
  // 2 hits of XXc at second pos; nothing at first pos (since that requires three AAA)
  testCase(t, "XXXc", "adc@1,  dac@1", needles, __LINE__);
  
  ///
  ///  new peptide DB
  ///
  
  t = ACTrie(2, 0);
  needles = {"eq", "nd", "llll"};
  t.addNeedlesAndCompress(needles);
  
  ///
  /// hits across the protein
  ///
  // B = D|N,  Z = E|Q
  // both match XX@1, eq matches ZZ, nd matches BB
  testCase(t, "aXXaBBkkZZlllllk", "nd@1, nd@4, eq@1, eq@8, llll@10, llll@11   ", needles, __LINE__);
  
  ///
  /// mismatches
  ///
  ///
  /// same, but with mm's allowed (but not sufficient)   
  ///
  t = ACTrie(0, 1);
  needles = {"acd", "adc", "cad", "cda", "dac", "dca"};
  t.addNeedlesAndCompress(needles);
  testCase(t, "aaaIIcccIIddd", "", needles, __LINE__);
  
  ///
  /// full usage of mm's
  ///
  t = ACTrie(0, 3);
  t.addNeedlesAndCompress(needles);
  testCase(t, "mmmm",
              "  dac@0,  cad@0,  cda@0,  dca@0,  adc@0,  acd@0"  // all six hits, found at first position
              ", dac@1,  cad@1,  cda@1,  dca@1,  adc@1,  acd@1"  // all six hits, found at second position
              , needles, __LINE__);

  ///
  /// with prefix
  ///
  t.setMaxMMCount(2);

  testCase(t, "aMMM",
           "acd@0,  adc@0", // 2 hits of aXX at first pos
           needles, __LINE__);

  ///
  /// with prefix and B 
  ///
  t.setMaxAAACount(1);

  testCase(t, "aMMB",
           "  adc@0,  acd@0"  // 2 hits of aXx at first pos
           ", cad@1,  acd@1"  // 2 hits of XXB, found at second position
           , needles, __LINE__);

  ///
  ///  new peptide DB
  ///
  t = ACTrie(1, 1);
  needles = {"eq", "nd", "llll"};
  t.addNeedlesAndCompress(needles);

  ///
  /// hits across the protein
  ///
  testCase(t, "aXXaBBkkZZlllllk",
           "nd@0, nd@1, nd@2, nd@3, nd@4, nd@5, eq@0, eq@1, eq@2, eq@7, eq@8, eq@9, llll@9, llll@10, llll@11, llll@12    " 
            //   nd matches all positions of 'aXXaBk';;  eq matches 'aXXa' and 'kZZl' ;; llll matches 'Zlllllk' 
           ,
           needles, __LINE__);


  ///
  /// matching Peptides WITH AAA's in them (should just be matched without digging into AAA/MM reserves)
  ///   
  t = ACTrie(0, 0);
  needles = {"acb", "abc", "cda", "bac", "anc", "acn", "dad"};
  t.addNeedlesAndCompress(needles);
  testCase(t, "baxyacbIIabcIIbac", "acb@4, abc@9, bac@14", needles, __LINE__);

  t = ACTrie(1, 0);
  // B = D|N,  Z = E|Q
  t.addNeedlesAndCompress(needles);
  testCase(t, "baxyacbIIabcIIbac", "acb@4, abc@9, bac@0, bac@14, anc@9, acn@4", needles, __LINE__);


  t = ACTrie(2, 0);
  // B = D|N,  Z = E|Q
  needles = {"dad", "bax", "bac", "anc", "acn"};
  t.addNeedlesAndCompress(needles);
  testCase(t, "baxyacbIIabcIIbac", "dad@0, bax@0, bac@0, bac@14, anc@9, acn@4", needles, __LINE__);

  t = ACTrie(2, 2);
  // B = D|N,  Z = E|Q
  needles = {"dady", "baxy", "iibac", "ancii", "yaknif"};
  t.addNeedlesAndCompress(needles);
  testCase(t, "baxyacbIIabcIIbac", "dady@0, baxy@0, dady@8, ancii@4, ancii@9, iibac@1, iibac@7, iibac@12, yaknif@3", needles, __LINE__);

  t = ACTrie(0, 0);
  needles = {"PEPTIDER", "XXXBEBEAR"};
  t.addNeedlesAndCompress(needles);
  testCase(t, "PEPTIDERXXXBEBEAR", "PEPTIDER@0, XXXBEBEAR@8", needles, __LINE__);
  
  ///
  /// TEST if offsets into proteins are correct in the presence of non-AA characters like '*'
  ///        NOTE: offsets will be incorrect if a hit overlaps with a '*', since the trie only knows the length of a hit and the end position
  ///              in the protein, thus computing the start will be off by the amount of '*'s
  t = ACTrie(0, 0);
  needles = {"MLTEAEK"};
  t.addNeedlesAndCompress(needles);
  testCase(t, "*MLTEAXK*", "", needles, __LINE__);

  t = ACTrie(1, 0);
  needles = {"MLTEAEK"};
  t.addNeedlesAndCompress(needles);
  testCase(t, "*MLTEAXK*", "MLTEAEK@1", needles, __LINE__);

  ///
  /// test if spawn does not report hits which do not cover its first AAA
  /// 
  t = ACTrie(4, 0);
  needles = {"MDDDEADC", "MDD", "DD", "DEADC"};
  t.addNeedlesAndCompress(needles);
  testCase(t, "MBBDEABCRAFG", "MDDDEADC@0, MDD@0, DD@1, DD@2, DEADC@3", needles, __LINE__);
             //MDDDEADC
}
END_SECTION

START_SECTION(void getAllHits(ACTrieState& state) const)
  NOT_TESTABLE // tested above
END_SECTION

/////////////////////////////////////
//// testing ACTrieState
/////////////////////////////////////

START_SECTION(void setQuery(const std::string& haystack))
  ACTrieState state;
  string q("PFAINGER");
  state.setQuery(q);
  TEST_EQUAL(state.getQuery(), q);
  TEST_EQUAL(state.hits.size(), 0);
  TEST_EQUAL(state.tree_pos(), 0);
  TEST_EQUAL(state.textPosIt() == state.getQuery().cbegin(), true);
  TEST_EQUAL(state.spawns.empty(), true);
  TEST_EQUAL(state.textPos(), 0);
END_SECTION

START_SECTION(size_t textPos() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(std::string::const_iterator textPosIt() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(const std::string& getQuery() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(AA nextValidAA())
  ACTrieState state;
  {
    constexpr string_view sv = "PFBNX";
    state.setQuery(string(sv.data()));
    checkAAIterator<5>({AA('P'), AA('F'), AA('B'), AA('N'), AA('X')}, {1, 2, 3, 4, 5}, state);
  }
  {
    constexpr string_view sv = "?X*B**";
    state.setQuery(string(sv.data()));
    TEST_EQUAL(state.getQuery(), sv);
    checkAAIterator<2>({AA('X'), AA('B')}, {2, 4}, state);
  }
  {
    state.setQuery("");
    TEST_EQUAL(state.nextValidAA().isValid(), false);
  }
END_SECTION 

START_SECTION(AA nextValidAA())
  NOT_TESTABLE // tested above
END_SECTION
  

/////////////////////////////////////
//// testing AA
/////////////////////////////////////

START_SECTION(constexpr AA())
{
  // make sure ctor is constexpr
  static_assert(AA('?').isValid() == false);

  static_assert(AA('?')() == CharToAA[(unsigned char)'?']);

  static_assert(AA('G') <= AA('B'));
  
  static_assert(AA('B').isAmbiguous());
  static_assert(AA('J').isAmbiguous());
  static_assert(AA('Z').isAmbiguous());
  static_assert(AA('X').isAmbiguous());
  static_assert(AA('$').isAmbiguous());

  static_assert(AA('B')++ == AA('B'));
  static_assert(++AA('B') == AA('J'));

  static_assert((AA('B') - AA('B'))() == 0);
  static_assert((AA('J') - AA('B'))() == 1);
  static_assert((AA('Z') - AA('B'))() == 2);
  static_assert((AA('X') - AA('B'))() == 3);

  for (char c = 'A'; c <= 'Z'; ++c)
    assert(AA(c).isValidForPeptide());

  static_assert(!AA('?').isValidForPeptide());
  static_assert(!AA('$').isValidForPeptide());
  static_assert(!AA(' ').isValidForPeptide());
  static_assert(!AA('*').isValidForPeptide());
  static_assert(!AA('3').isValidForPeptide());
  static_assert(!AA('#').isValidForPeptide());
}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
