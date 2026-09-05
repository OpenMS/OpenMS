// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/MS1LabelState.h>
#include <OpenMS/METADATA/PeptideHit.h>

using namespace OpenMS;

START_TEST(MS1LabelState, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((const std::vector<std::string>& keys()))
{
  const auto& keys = MS1LabelState::keys();
  TEST_EQUAL(keys.size(), 3)
  TEST_EQUAL(keys[0], MS1LabelState::LABELED_SEQUENCE)
  TEST_EQUAL(keys[1], MS1LabelState::REMOVED_LABELS)
  TEST_EQUAL(keys[2], MS1LabelState::CHANNEL)
  TEST_EQUAL(MS1LabelState::LABELED_SEQUENCE, "MS1Label:labeled_sequence")
  TEST_EQUAL(MS1LabelState::REMOVED_LABELS, "MS1Label:removed_labels")
  TEST_EQUAL(MS1LabelState::CHANNEL, "MS1Label:channel")
}
END_SECTION

START_SECTION((bool hasMatchedSequence(const PeptideHit& hit)))
{
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  TEST_FALSE(MS1LabelState::hasMatchedSequence(hit))
  hit.setMetaValue(MS1LabelState::LABELED_SEQUENCE, "PEPTIDEK(Label:13C(6)15N(2))");
  TEST_TRUE(MS1LabelState::hasMatchedSequence(hit))
}
END_SECTION

START_SECTION((AASequence matchedSequence(const PeptideHit& hit)))
{
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  // no label state recorded: the sequence itself
  TEST_EQUAL(MS1LabelState::matchedSequence(hit).toString(), "PEPTIDEK")

  // reduced to the peptide identity: the peptidoform as matched
  hit.setMetaValue(MS1LabelState::LABELED_SEQUENCE, "PEPTIDEK(Label:13C(6)15N(2))");
  const AASequence matched = MS1LabelState::matchedSequence(hit);
  TEST_EQUAL(matched.toString(), "PEPTIDEK(Label:13C(6)15N(2))")
  TEST_TRUE(matched.isModified())
  TEST_EQUAL(matched.toUnmodifiedString(), "PEPTIDEK")
  // the label mass is part of the matched form, not of the identity
  TEST_REAL_SIMILAR(matched.getMonoWeight() - hit.getSequence().getMonoWeight(), 8.0141988132)
}
END_SECTION

START_SECTION((PeptideHit withMatchedSequence(const PeptideHit& hit)))
{
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(2);
  hit.setScore(0.01);

  // without a label state the copy is unchanged
  PeptideHit same = MS1LabelState::withMatchedSequence(hit);
  TEST_EQUAL(same.getSequence().toString(), "PEPTIDEK")
  TEST_EQUAL(same.getCharge(), 2)

  hit.setMetaValue(MS1LabelState::LABELED_SEQUENCE, "PEPTIDEK(Label:13C(6)15N(2))");
  hit.setMetaValue(MS1LabelState::REMOVED_LABELS, "Lys8");
  hit.setMetaValue(MS1LabelState::CHANNEL, 2);
  PeptideHit matched = MS1LabelState::withMatchedSequence(hit);
  TEST_EQUAL(matched.getSequence().toString(), "PEPTIDEK(Label:13C(6)15N(2))")
  // everything else, the label state included, is carried over; the original is untouched
  TEST_EQUAL(matched.getCharge(), 2)
  TEST_REAL_SIMILAR(matched.getScore(), 0.01)
  TEST_EQUAL(matched.getMetaValue(MS1LabelState::REMOVED_LABELS).toString(), "Lys8")
  TEST_EQUAL((int)matched.getMetaValue(MS1LabelState::CHANNEL), 2)
  TEST_EQUAL(hit.getSequence().toString(), "PEPTIDEK")
}
END_SECTION

START_SECTION((bool hasMatchedSequence(const PeptideHit& hit, const Keys& keys)))
{
  // The index overloads are what every exporter uses, so they are exercised on their own here.
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setMetaValue(MS1LabelState::LABELED_SEQUENCE, "PEPTIDEK(Label:13C(6)15N(2))");

  // Resolved after the keys are registered by the setMetaValue above, as the docs require.
  const MS1LabelState::Keys keys;
  TEST_TRUE(MS1LabelState::hasMatchedSequence(hit, keys))
  TEST_EQUAL(MS1LabelState::matchedSequence(hit, keys).toString(), "PEPTIDEK(Label:13C(6)15N(2))")

  // A hit without the key, read through the same resolved indices
  PeptideHit plain;
  plain.setSequence(AASequence::fromString("PEPTIDER"));
  TEST_FALSE(MS1LabelState::hasMatchedSequence(plain, keys))
  TEST_EQUAL(MS1LabelState::matchedSequence(plain, keys).toString(), "PEPTIDER")

  // The "no hit in this process has used the key yet" sentinel: an all-UInt(-1) Keys must read as
  // absent rather than indexing the registry with it.
  MS1LabelState::Keys unregistered = keys;
  unregistered.labeled_sequence = static_cast<UInt>(-1);
  TEST_FALSE(MS1LabelState::hasMatchedSequence(hit, unregistered))
  TEST_EQUAL(MS1LabelState::matchedSequence(hit, unregistered).toString(), "PEPTIDEK")
}
END_SECTION

START_SECTION(([EXTRA] an unusable labeled_sequence falls back to the sequence of the hit))
{
  // The value survives consensusXML/idXML as a plain UserParam, so it can arrive edited, foreign or
  // stale. None of these may throw: the exporters stream, so a throw abandons a half-written file.
  const AASequence own = AASequence::fromString("PEPTIDEK");

  // empty: parses to an empty AASequence rather than throwing, so it has to be rejected up front or
  // it silently replaces the sequence with nothing and the row is dropped
  PeptideHit empty_value;
  empty_value.setSequence(own);
  empty_value.setMetaValue(MS1LabelState::LABELED_SEQUENCE, "");
  TEST_FALSE(MS1LabelState::hasMatchedSequence(empty_value))
  TEST_EQUAL(MS1LabelState::matchedSequence(empty_value).toString(), "PEPTIDEK")
  TEST_EQUAL(MS1LabelState::withMatchedSequence(empty_value).getSequence().toString(), "PEPTIDEK")

  // unknown modification in round brackets -> Exception::InvalidValue
  PeptideHit unknown_mod;
  unknown_mod.setSequence(own);
  unknown_mod.setMetaValue(MS1LabelState::LABELED_SEQUENCE, "PEPTIDEK(NoSuchModification)");
  TEST_TRUE(MS1LabelState::hasMatchedSequence(unknown_mod))
  TEST_EQUAL(MS1LabelState::matchedSequence(unknown_mod).toString(), "PEPTIDEK")

  // square brackets -> Exception::ConversionError, not ParseError or InvalidValue; this is the form
  // the QPX peptidoform column emits, so it is the likeliest bad value to arrive here
  PeptideHit bracket_mod;
  bracket_mod.setSequence(own);
  bracket_mod.setMetaValue(MS1LabelState::LABELED_SEQUENCE, "PEPTIDEK[NotANumber]");
  TEST_TRUE(MS1LabelState::hasMatchedSequence(bracket_mod))
  TEST_EQUAL(MS1LabelState::matchedSequence(bracket_mod).toString(), "PEPTIDEK")

  // not a sequence at all -> Exception::ParseError
  PeptideHit junk;
  junk.setSequence(own);
  junk.setMetaValue(MS1LabelState::LABELED_SEQUENCE, "not a peptide!");
  TEST_EQUAL(MS1LabelState::matchedSequence(junk).toString(), "PEPTIDEK")

  // the index overload falls back the same way
  const MS1LabelState::Keys keys;
  TEST_EQUAL(MS1LabelState::matchedSequence(bracket_mod, keys).toString(), "PEPTIDEK")
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
