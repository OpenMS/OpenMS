// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include "ModifiedNASequenceGenerator.h"
#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
///////////////////////////

#include <string>
#include <set>
#include <cmath>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

using namespace OpenMS;
using namespace std;

START_TEST(ModifiedNASequenceGenerator, "$Id$")

const RibonucleotideDB* db = RibonucleotideDB::getInstance();

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

  // if a residue is already modified, it must not be changed again;
  // this is relevant for NASE where pre-digest (cleavage-sensitive) mods
  // consume part of the variable-mod budget.
  ams.clear();
  var_mods.clear();
  var_mods.insert(db->getRibonucleotide("m3U"));
  var_mods.insert(db->getRibonucleotide("s4U"));
  sequence = NASequence::fromString("A[m2,2G]UA");

  ModifiedNASequenceGenerator::applyVariableModifications(var_mods, sequence, 0, ams, true);
  TEST_EQUAL(ams.size(), 1);
  TEST_STRING_EQUAL(ams[0].toString(), sequence.toString());

  ams.clear();
  ModifiedNASequenceGenerator::applyVariableModifications(var_mods, sequence, 1, ams, true);
  // original + two one-site U modifications at the unmodified U position
  TEST_EQUAL(ams.size(), 3);
  set<string> seqs;
  for (const NASequence& s : ams) seqs.insert(s.toString());
  TEST_EQUAL(seqs.count(sequence.toString()), 1);
  TEST_EQUAL(seqs.count(NASequence::fromString("A[m2,2G][m3U]A").toString()), 1);
  TEST_EQUAL(seqs.count(NASequence::fromString("A[m2,2G][s4U]A").toString()), 1);
}
END_SECTION


START_SECTION(([EXTRA] applyVariableModifications includes terminal variants in the one-modification path))
{
  const auto* db = RibonucleotideDB::getInstance();
  std::set<const Ribonucleotide*> mods{db->getRibonucleotide("5'-p"), db->getRibonucleotide("3'-p")};
  std::vector<NASequence> variants;
  ModifiedNASequenceGenerator::applyVariableModifications(mods, NASequence::fromString("ACGU"), 1, variants);
  std::set<std::string> sequences;
  for (const auto& variant : variants) sequences.insert(variant.toString());
  TEST_EQUAL(sequences.size(), 3)
  TEST_EQUAL(sequences.count("pACGU"), 1)
  TEST_EQUAL(sequences.count("ACGUp"), 1)
  variants.clear();
  ModifiedNASequenceGenerator::applyVariableModifications({}, NASequence(), 1, variants);
  TEST_EQUAL(variants.size(), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


START_SECTION(([EXTRA] constrained placement preserves exhaustive compatibility, budgets and terminals))
{
  const std::set<const Ribonucleotide*> mods {db->getRibonucleotide("mA"),  db->getRibonucleotide("Am"),    db->getRibonucleotide("m3U"),
                                              db->getRibonucleotide("s4U"), db->getRibonucleotide("m2,2G"), db->getRibonucleotide("5'-p"),
                                              db->getRibonucleotide("3'-p")};
  for (const auto& text : {"", "ACGU", "pA[mA]GUp", "AC[m2,2G]UA"})
  {
    const auto source = NASequence::fromString(text);
    for (Size budget : {0, 1, 2, 3})
    {
      for (bool keep : {false, true})
      {
        std::vector<NASequence> reference, actual;
        ModifiedNASequenceGenerator::applyVariableModifications(mods, source, budget, reference, keep);
        ModifiedNASequenceGenerator::PlacementConstraints constraints;
        constraints.accepts = [](double, double, const auto&) { return true; };
        const auto statistics = ModifiedNASequenceGenerator::applyVariableModificationsConstrained(mods, source, budget, actual, constraints, keep);
        std::sort(reference.begin(), reference.end());
        std::sort(actual.begin(), actual.end());
        TEST_TRUE(actual == reference)
        TEST_FALSE(statistics.enumeration_fallback)
        TEST_EQUAL(statistics.emitted_variants, actual.size())
      }
    }
  }
}
END_SECTION

START_SECTION(([EXTRA] conservative mass bounds prune partial placements and retain negative mass shifts))
{
  Ribonucleotide loss = *db->getRibonucleotide("A");
  loss.setCode("test_loss");
  loss.setOrigin('A');
  loss.setFormula(loss.getFormula() - EmpiricalFormula("H2"));
  const std::set<const Ribonucleotide*> mods {&loss,
                                              db->getRibonucleotide("mA"),
                                              db->getRibonucleotide("Am"),
                                              db->getRibonucleotide("m3U"),
                                              db->getRibonucleotide("5'-p"),
                                              db->getRibonucleotide("3'-p")};
  const auto source = NASequence::fromString("AAAAUU");
  std::vector<NASequence> reference;
  ModifiedNASequenceGenerator::applyVariableModifications(mods, source, 3, reference);
  for (bool average : {false, true})
  {
    const auto massOf = [average](const NASequence& sequence) { return average ? sequence.getAverageWeight() : sequence.getMonoWeight(); };
    NASequence desired = source;
    desired.set(1, &loss);
    const double target = massOf(desired);
    ModifiedNASequenceGenerator::PlacementConstraints constraints;
    constraints.use_average_mass = average;
    constraints.accepts = [target](double low, double high, const auto&) { return low <= target && target <= high; };
    std::vector<NASequence> actual, expected;
    const auto statistics = ModifiedNASequenceGenerator::applyVariableModificationsConstrained(mods, source, 3, actual, constraints);
    for (const auto& candidate : reference)
    {
      if (std::abs(massOf(candidate) - target) < 1e-7) { expected.push_back(candidate); }
    }
    std::sort(actual.begin(), actual.end());
    std::sort(expected.begin(), expected.end());
    TEST_TRUE(actual == expected)
    TEST_TRUE(statistics.pruned_branches > 0)
    TEST_TRUE(statistics.emitted_variants < reference.size())
    TEST_FALSE(statistics.enumeration_fallback)
  }
  ModifiedNASequenceGenerator::PlacementConstraints bounded;
  bounded.max_nodes = 1;
  bounded.accepts = [](double, double, const auto&) { return true; };
  std::vector<NASequence> output {NASequence::fromString("C")};
  const auto statistics = ModifiedNASequenceGenerator::applyVariableModificationsConstrained(mods, source, 3, output, bounded);
  TEST_TRUE(statistics.enumeration_fallback)
  TEST_EQUAL(output.size(), reference.size() + 1)
  TEST_STRING_EQUAL(output.front().toString(), "C")
  output.erase(output.begin());
  std::sort(output.begin(), output.end());
  std::sort(reference.begin(), reference.end());
  TEST_TRUE(output == reference)
}
END_SECTION

END_TEST


