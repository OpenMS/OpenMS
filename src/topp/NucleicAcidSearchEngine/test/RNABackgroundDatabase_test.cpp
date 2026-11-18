// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#include "RNABackgroundDatabase.h"

#include <OpenMS/CHEMISTRY/RNaseDB.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CONCEPT/ClassTest.h>

#include <algorithm>
#include <set>

using namespace OpenMS;

namespace
{
const DigestionEnzymeRNA& t1()
{
  return *RNaseDB::getInstance()->getEnzyme("RNase_T1");
}

std::multiset<char> composition(const std::string& s)
{
  return std::multiset<char>(s.begin(), s.end());
}

const auto nothing = [](const std::string&) { return false; };
} // namespace

START_TEST(RNABackgroundDatabase, "$Id$")

START_SECTION(isFixed: cleavage residues of the enzyme)
{
  RNABackgroundDatabase db(t1(), 5, 1);
  auto* ribos = RibonucleotideDB::getInstance();
  TEST_EQUAL(db.isFixed(*ribos->getRibonucleotide("G")), true)
  TEST_EQUAL(db.isFixed(*ribos->getRibonucleotide("A")), false)
  TEST_EQUAL(db.isFixed(*ribos->getRibonucleotide("C")), false)
  TEST_EQUAL(db.isFixed(*ribos->getRibonucleotide("U")), false)
  // RNase T1's 'cuts before' pattern ('.*') matches every nucleotide and fixes none
}
END_SECTION

START_SECTION(productKeys: products without missed cleavages and at least min_length residues)
{
  RNABackgroundDatabase db(t1(), 5, 1);
  std::vector<std::string> keys = db.productKeys("pGCAUCCAUGAGUUCCACCUGCCA");
  // products: G, CAUCCAUG, AG, UUCCACCUG, CCA
  TEST_EQUAL(keys.size(), 2)
  TEST_STRING_EQUAL(keys[0], "CAUCCAUG")
  TEST_STRING_EQUAL(keys[1], "UUCCACCUG")
}
END_SECTION

START_SECTION(shuffleProduct: same composition with cleavage residues and the last residue in place)
{
  RNABackgroundDatabase db(t1(), 5, 7);
  NASequence product = NASequence::fromString("UUCCAGCACCUG");
  auto shuffled = db.shuffleProduct(product, nothing);
  TEST_EQUAL(bool(shuffled), true)
  const std::string s = RNABackgroundDatabase::key(*shuffled);
  TEST_EQUAL(s.size(), 12)
  TEST_EQUAL(s != "UUCCAGCACCUG", true)
  TEST_EQUAL(composition(s) == composition("UUCCAGCACCUG"), true)
  TEST_EQUAL(s[5], 'G')
  TEST_EQUAL(s[11], 'G')
  // deterministic per seed and product
  RNABackgroundDatabase same(t1(), 5, 7);
  TEST_STRING_EQUAL(RNABackgroundDatabase::key(*same.shuffleProduct(product, nothing)), s)
  // forbidden shuffles are avoided
  auto other = db.shuffleProduct(product, [&](const std::string& k) { return k == s; });
  TEST_EQUAL(bool(other), true)
  TEST_EQUAL(RNABackgroundDatabase::key(*other) != s, true)
  // no shuffle if fewer than two movable residues differ in position
  TEST_EQUAL(bool(db.shuffleProduct(NASequence::fromString("AAAAG"), nothing)), false)
  // ... or if every shuffle is forbidden
  TEST_EQUAL(bool(db.shuffleProduct(NASequence::fromString("ACAAG"), [](const std::string&) { return true; })), false)
}
END_SECTION

START_SECTION(rewrite: selected products are replaced consistently across sequences)
{
  RNABackgroundDatabase db(t1(), 5, 3);
  auto select = [](const std::string& k) { return k == "UUCCACCUG"; };
  const std::string a = db.rewrite("pGCAUCCAUGAGUUCCACCUGCCA", select, nothing);
  const std::string b = db.rewrite("pGUUCCACCUGCCA", select, nothing);
  TEST_EQUAL(a.substr(0, 12), "pGCAUCCAUGAG") // unselected products unchanged, 5' phosphate kept
  TEST_EQUAL(a.size(), std::string("pGCAUCCAUGAGUUCCACCUGCCA").size())
  const std::string replaced = a.substr(12, 9);
  TEST_EQUAL(replaced != "UUCCACCUG", true)
  TEST_EQUAL(replaced.back(), 'G')
  TEST_STRING_EQUAL(b, "pG" + replaced + "CCA") // the same product gets the same replacement
  TEST_EQUAL(db.rewritten(), 1)
  TEST_EQUAL(db.failed(), 0)
  TEST_EQUAL(db.replacements().count(replaced), 1)
  // a later replacement never repeats an earlier one
  auto all = [](const std::string&) { return true; };
  RNABackgroundDatabase db2(t1(), 5, 3);
  const std::string c = db2.rewrite("UUCCACCUG", all, nothing);
  const std::string d = db2.rewrite("CUUCACCUG", all, [&](const std::string& k) { return k == "CUUCACCUG"; });
  TEST_EQUAL(c != d, true)
  // short products are never selected
  TEST_STRING_EQUAL(db2.rewrite("pGCAG", all, nothing), "pGCAG")
  // ... unless 'shuffle_short' is set: then they are shuffled too (here 'pAUCG')
  const std::string short_shuffled = db2.rewrite("pAUCGCUCCG", all, nothing, true);
  TEST_EQUAL(short_shuffled != "pAUCGCUCCG", true)
  TEST_EQUAL(short_shuffled.size(), 10)
  TEST_EQUAL(short_shuffled[4], 'G')
  TEST_EQUAL(composition(short_shuffled) == composition("pAUCGCUCCG"), true)
}
END_SECTION

START_SECTION(rewrite: missed-cleavage combinations of short products are checked)
{
  // review case (OpenMS PR #10941): all products of these targets are shorter than 5, but with one missed cleavage
  // each target is a searchable 6-mer; a naive shuffle maps each target to the other one
  const std::vector<std::string> targets = {"ACGUAG", "CAGAUG"};
  RNABackgroundDatabase db(t1(), 5, 1);
  std::set<std::string> target_keys;
  for (const auto& target : targets)
  {
    for (const auto& k : db.productKeys(target, 1)) target_keys.insert(k);
  }
  TEST_EQUAL(target_keys.size(), 2) // ACGUAG, CAGAUG
  auto is_target = [&](const std::string& k) { return target_keys.count(k) > 0; };
  auto all = [](const std::string&) { return true; };
  const std::vector<std::string> decoys = db.rewrite(targets, all, is_target, true, 1);
  TEST_EQUAL(decoys.size(), 2)
  for (const auto& decoy : decoys)
  {
    for (const auto& k : db.productKeys(decoy, 1)) TEST_EQUAL(is_target(k), false)
    TEST_EQUAL(decoy[2], 'G') // cleavage residues stay
    TEST_EQUAL(decoy[5], 'G')
  }
  TEST_EQUAL(db.failedCombinations(), 0)
  // deterministic
  RNABackgroundDatabase same(t1(), 5, 1);
  const std::vector<std::string> again = same.rewrite(targets, all, is_target, true, 1);
  TEST_EQUAL(again == decoys, true)
}
END_SECTION

START_SECTION(productKeys: missed-cleavage combinations)
{
  RNABackgroundDatabase db(t1(), 5, 1);
  // products: ACG, UAG, CCCCCG, A
  std::vector<std::string> keys = db.productKeys("ACGUAGCCCCCGA", 1);
  std::set<std::string> expected = {"ACGUAG", "UAGCCCCCG", "CCCCCG", "CCCCCGA"};
  TEST_EQUAL(std::set<std::string>(keys.begin(), keys.end()) == expected, true)
}
END_SECTION

END_TEST
