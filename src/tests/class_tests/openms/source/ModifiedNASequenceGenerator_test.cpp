// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ModifiedNASequenceGenerator.h>
#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
///////////////////////////

#include <string>

using namespace OpenMS;
using namespace std;

START_TEST(ModifiedNASequenceGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
RibonucleotideDB* db = RibonucleotideDB::getInstance();

START_SECTION((static void applyFixedModifications(const std::set<ModifiedNASequenceGenerator::ConstRibonucleotidePtr>& fixed_mods, NASequence& sequence)))
{
  set<ModifiedNASequenceGenerator::ConstRibonucleotidePtr> fixed_mods;

  // query modified ribos by code
  vector<string> fixed_mods_code = {"s4U"};  // 4-thiouridine
  for (auto const & f : fixed_mods_code) { fixed_mods.insert(db->getRibonucleotide(f)); }

  NASequence sequence = NASequence::fromString("AUAUAUA");

  ModifiedNASequenceGenerator::applyFixedModifications(fixed_mods, sequence);

  TEST_STRING_EQUAL(sequence.toString(), "A[s4U]A[s4U]A[s4U]A");

  // additional check if internal representation equal
  NASequence sequence2 = NASequence::fromString("A[s4U]A[s4U]A[s4U]A");
  TEST_EQUAL(sequence, sequence2);
}
END_SECTION

START_SECTION(static void applyVariableModifications(const std::set<ConstRibonucleotidePtr>& var_mods, const NASequence& seq, Size max_variable_mods_per_NASequence, std::vector<NASequence>& all_modified_NASequences, bool keep_original = true))
{
  set<ModifiedNASequenceGenerator::ConstRibonucleotidePtr> var_mods;
  // query modified ribos by code
  vector<string> mods_code = {"m3U", "s4U"};  // 3-methyluridine, 4-thiouridine
  for (auto const & f : mods_code) { var_mods.insert(db->getRibonucleotide(f)); }

  NASequence sequence = NASequence::fromString("AUAUAUA");
  vector<NASequence> ams;

  // (1) Add at most one modification. (true) return the unmodified version
  ModifiedNASequenceGenerator::applyVariableModifications(var_mods, sequence, 1, ams, true);

  TEST_EQUAL(ams.size(), 7);
  TEST_STRING_EQUAL(ams[0].toString(), NASequence::fromString("AUAUAUA").toString());
  // the order of "m3U" and "s4U" in "var_mods" is unclear (pointers ordered by
  // address) and determines the order of the result ("ams") - need to sort:
  sort(++ams.begin(), ams.end());
  TEST_STRING_EQUAL(ams[1].toString(), NASequence::fromString("AUAUA[m3U]A").toString());
  TEST_STRING_EQUAL(ams[2].toString(), NASequence::fromString("AUAUA[s4U]A").toString());
  TEST_STRING_EQUAL(ams[3].toString(), NASequence::fromString("AUA[m3U]AUA").toString());
  TEST_STRING_EQUAL(ams[4].toString(), NASequence::fromString("AUA[s4U]AUA").toString());
  TEST_STRING_EQUAL(ams[5].toString(), NASequence::fromString("A[m3U]AUAUA").toString());
  TEST_STRING_EQUAL(ams[6].toString(), NASequence::fromString("A[s4U]AUAUA").toString());

  ams.clear();
  // (1) Add at most one modification. (false) without the unmodified version
  ModifiedNASequenceGenerator::applyVariableModifications(var_mods, sequence, 1, ams, false);

  TEST_EQUAL(ams.size(), 6); // same as before but now without the unmodified version

  ams.clear();
  // (3) Add at most three modification. (true) with the unmodified version
  ModifiedNASequenceGenerator::applyVariableModifications(var_mods, sequence, 3, ams, true);
  TEST_EQUAL(ams.size(), 3*3*3); // 3^3 sequences expected

  // test modification of A and U
  ams.clear();
  var_mods.clear();
  mods_code = {"s4U", "m3U", "m1A"};
  for (auto const & f : mods_code) { var_mods.insert(db->getRibonucleotide(f)); }

  ModifiedNASequenceGenerator::applyVariableModifications(var_mods, sequence, 7, ams, true);
  TEST_EQUAL(ams.size(), 3*3*3*2*2*2*2); // 3^3 combinations for U times 2^4 for A
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



