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
  TEST_EQUAL(keys[2], MS1LabelState::LABEL_CHANNEL)
  TEST_EQUAL(MS1LabelState::LABELED_SEQUENCE, "labeled_sequence")
  TEST_EQUAL(MS1LabelState::REMOVED_LABELS, "removed_labels")
  TEST_EQUAL(MS1LabelState::LABEL_CHANNEL, "label_channel")
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
  hit.setMetaValue(MS1LabelState::LABEL_CHANNEL, 2);
  PeptideHit matched = MS1LabelState::withMatchedSequence(hit);
  TEST_EQUAL(matched.getSequence().toString(), "PEPTIDEK(Label:13C(6)15N(2))")
  // everything else, the label state included, is carried over; the original is untouched
  TEST_EQUAL(matched.getCharge(), 2)
  TEST_REAL_SIMILAR(matched.getScore(), 0.01)
  TEST_EQUAL(matched.getMetaValue(MS1LabelState::REMOVED_LABELS).toString(), "Lys8")
  TEST_EQUAL((int)matched.getMetaValue(MS1LabelState::LABEL_CHANNEL), 2)
  TEST_EQUAL(hit.getSequence().toString(), "PEPTIDEK")
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
