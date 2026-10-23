// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>

using namespace OpenMS;
using Adapter = IdentificationDataAdapter;

namespace
{
ProteinIdentification protein()
{
  ProteinIdentification result;
  result.setIdentifier("search");
  return result;
}
PeptideIdentification peptide()
{
  PeptideIdentification result;
  result.setIdentifier("search");
  result.setScoreType("score");
  result.setHits({PeptideHit(3, 1, 2, AASequence::fromString("PEPTIDE")), PeptideHit(2, 2, 2, AASequence::fromString("OTHER"))});
  return result;
}
} // namespace

START_TEST(IdentificationDataWorkflow, "$Id$")

START_SECTION([EXTRA] feature and subordinate links are live while measured values survive filtering)
{
  FeatureMap map;
  map.setProteinIdentifications({protein()});
  Feature feature;
  feature.setIntensity(1200.5);
  feature.setRT(25);
  feature.getPeptideIdentifications().push_back(peptide());
  Feature subordinate;
  subordinate.setIntensity(42);
  subordinate.getPeptideIdentifications().push_back(peptide());
  feature.getSubordinates().push_back(subordinate);
  map.push_back(feature);
  map.getUnassignedPeptideIdentifications().push_back(peptide());
  auto imported = Adapter::fromFeatureMap(map);
  TEST_EQUAL(imported.associations.size(), 3)
  TEST_EQUAL(imported.associations[1].feature_path.size(), 2)
  TEST_TRUE(imported.associations[2].unassigned)
  auto& run = imported.data.getRun("search");
  run.filterMatches([](const auto& match) { return match.representation == "PEPTIDE"; });
  auto links = imported.associations;
  TEST_EXCEPTION(Exception::InvalidParameter, Adapter::reconcileAssociations(imported.data, links, Adapter::MissingLinkPolicy::REJECT))
  TEST_EQUAL(links[0].matches.size(), 2)
  TEST_EQUAL(Adapter::reconcileAssociations(imported.data, links, Adapter::MissingLinkPolicy::PRUNE), 3)
  Adapter::ExportOptions options;
  auto losses = Adapter::applyToFeatureMap(imported.data, links, map, options, Adapter::MissingLinkPolicy::REJECT);
  TEST_TRUE(losses.empty())
  TEST_REAL_SIMILAR(map[0].getIntensity(), 1200.5)
  TEST_REAL_SIMILAR(map[0].getSubordinates()[0].getIntensity(), 42)
  TEST_EQUAL(map[0].getPeptideIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(map[0].getSubordinates()[0].getPeptideIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(map.getUnassignedPeptideIdentifications()[0].getHits().size(), 1)
  run.eraseMatches([](const auto&) { return true; });
  Adapter::applyToFeatureMap(imported.data, links, map, options, Adapter::MissingLinkPolicy::PRUNE);
  TEST_EQUAL(map.size(), 1)
  TEST_EQUAL(map[0].getPeptideIdentifications().size(), 0)
  TEST_REAL_SIMILAR(map[0].getIntensity(), 1200.5)
}
END_SECTION

START_SECTION([EXTRA] consensus intensities and channels remain independent of identification links)
{
  ConsensusMap map;
  map.setProteinIdentifications({protein()});
  ConsensusFeature feature;
  feature.setIntensity(999);
  feature.getPeptideIdentifications().push_back(peptide());
  map.push_back(feature);
  map.getColumnHeaders()[0].filename = "TMT-channel";
  map.getColumnHeaders()[0].label = "126";
  auto imported = Adapter::fromConsensusMap(map);
  auto& run = imported.data.getRun("search");
  run.filterMatches([](const auto& match) { return match.representation == "PEPTIDE"; });
  Adapter::ExportOptions options;
  Adapter::applyToConsensusMap(imported.data, imported.associations, map, options, Adapter::MissingLinkPolicy::PRUNE);
  TEST_REAL_SIMILAR(map[0].getIntensity(), 999)
  TEST_EQUAL(map[0].getPeptideIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(map.getColumnHeaders()[0].label, "126")
  auto broken = imported.associations;
  broken[0].feature_path = {500};
  TEST_EXCEPTION(Exception::InvalidParameter, Adapter::applyToConsensusMap(imported.data, broken, map, options, Adapter::MissingLinkPolicy::PRUNE))
  TEST_REAL_SIMILAR(map[0].getIntensity(), 999)
  TEST_EQUAL(map[0].getPeptideIdentifications()[0].getHits().size(), 1)
}
END_SECTION

END_TEST
