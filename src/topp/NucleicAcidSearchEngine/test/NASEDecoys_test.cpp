// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#define OPENMS_NASE_CLASS_TEST
#include "../NucleicAcidSearchEngine.cpp"
#include <OpenMS/CONCEPT/ClassTest.h>

class TestNASEDecoys : public NucleicAcidSearchEngine
{
public:
  using NucleicAcidSearchEngine::reshuffleDecoyIfCollision_;
  using NucleicAcidSearchEngine::registerDigestedOligo_;
};

START_TEST(NASEDecoys, "$Id$")

START_SECTION(([EXTRA] collision retries are bounded and preserve chain termini))
{
  TestNASEDecoys tool;
  const NASequence original = NASequence::fromString("pACGUp");
  const set<NASequence> targets{original};
  bool resolved = false;
  Size attempts = 0;
  const auto first = tool.reshuffleDecoyIfCollision_(original, targets, 5, resolved, attempts);
  TEST_TRUE(resolved)
  TEST_EQUAL(attempts > 0 && attempts <= 5, true)
  TEST_EQUAL(targets.count(first), 0)
  TEST_EQUAL(first.getFivePrimeMod(), original.getFivePrimeMod())
  TEST_EQUAL(first.getThreePrimeMod(), original.getThreePrimeMod())
  const auto second = tool.reshuffleDecoyIfCollision_(original, targets, 5, resolved, attempts);
  TEST_EQUAL(first, second)
  tool.reshuffleDecoyIfCollision_(original, targets, 0, resolved, attempts);
  TEST_FALSE(resolved)
  TEST_EQUAL(attempts, 0)
  const auto unique = tool.reshuffleDecoyIfCollision_(first, targets, 0, resolved, attempts);
  TEST_TRUE(resolved)
  TEST_EQUAL(attempts, 0)
  TEST_EQUAL(unique, first)
  const auto homopolymer = NASequence::fromString("AAAA");
  const auto unresolved = tool.reshuffleDecoyIfCollision_(homopolymer, {homopolymer}, 5, resolved, attempts);
  TEST_FALSE(resolved)
  TEST_EQUAL(attempts, 5)
  TEST_EQUAL(unresolved, homopolymer)
}
END_SECTION

START_SECTION(([EXTRA] collision-free decoys retain annotated source coordinates))
{
  TestNASEDecoys tool;
  IdentificationData data;
  IdentificationData::ParentSequence target("TARGET", IdentificationData::MoleculeType::RNA, "pACGUp");
  IdentificationData::ParentSequence decoy("DECOY", IdentificationData::MoleculeType::RNA, "pACGUp");
  decoy.is_decoy = true;
  auto target_ref = data.registerParentSequence(target);
  auto decoy_ref = data.registerParentSequence(decoy);
  const auto parent = NASequence::fromString(target.sequence);
  set<NASequence> targets;
  Size collisions = 0, resolved = 0, unresolved = 0, exhausted = 0;
  tool.registerDigestedOligo_(data, target_ref, parent, {0, 4}, parent, targets, 5, collisions, resolved, unresolved, exhausted);
  tool.registerDigestedOligo_(data, decoy_ref, parent, {0, 4}, parent, targets, 5, collisions, resolved, unresolved, exhausted);
  TEST_EQUAL(collisions, 1)
  TEST_EQUAL(resolved, 1)
  TEST_EQUAL(unresolved, 0)
  TEST_EQUAL(data.getIdentifiedOligos().size(), 2)
  for (const auto& oligo : data.getIdentifiedOligos())
  {
    TEST_EQUAL(oligo.parent_matches.size(), 1)
    const auto& mapping = *oligo.parent_matches.begin();
    if (mapping.first->is_decoy)
    {
      TEST_NOT_EQUAL(oligo.sequence, parent)
      const auto& match = *mapping.second.begin();
      TEST_EQUAL(match.start_pos, 0)
      TEST_EQUAL(match.end_pos, 3)
      TEST_EQUAL(match.getMetaValue("reshuffled_decoy").toString(), "true")
    }
  }
  tool.registerDigestedOligo_(data, decoy_ref, parent, {0, 4}, parent, targets, 0, collisions, resolved, unresolved, exhausted);
  TEST_EQUAL(unresolved, 1)
  TEST_EQUAL(data.getIdentifiedOligos().size(), 2)
}
END_SECTION

END_TEST
