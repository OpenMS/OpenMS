// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>

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
  Feature feature;
  feature.setIntensity(1200.5);
  feature.setRT(25);
  Feature subordinate;
  subordinate.setIntensity(42);
  feature.getSubordinates().push_back(subordinate);
  map.push_back(feature);
  MapPeptideIdentifications ids;
  ids.proteins = {protein()};
  ids.features[{0}] = {peptide()};
  ids.features[{0, 0}] = {peptide()};
  ids.unassigned = {peptide()};
  IdentificationDataConverter::importMapIdentifications(map, ids);
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
  // filtering the matches of the map keeps the links of the remaining ones and the measured values
  map.eraseMatches([](const auto&, const auto&, const auto& match) { return match.representation != "PEPTIDE"; });
  const auto filtered = IdentificationDataConverter::exportMapIdentifications(map);
  TEST_REAL_SIMILAR(map[0].getIntensity(), 1200.5)
  TEST_REAL_SIMILAR(map[0].getSubordinates()[0].getIntensity(), 42)
  TEST_EQUAL(filtered.at({0})[0].getHits().size(), 1)
  TEST_EQUAL(filtered.at({0, 0})[0].getHits().size(), 1)
  TEST_EQUAL(filtered.unassigned[0].getHits().size(), 1)
  map.eraseIdentifications([](const auto&, const auto&) { return true; });
  TEST_EQUAL(map.size(), 1)
  TEST_EQUAL(IdentificationDataConverter::exportMapIdentifications(map).at({0}).size(), 0)
  TEST_REAL_SIMILAR(map[0].getIntensity(), 1200.5)
}
END_SECTION

START_SECTION([EXTRA] consensus intensities and channels remain independent of identification links)
{
  ConsensusMap map;
  ConsensusFeature feature;
  feature.setIntensity(999);
  map.push_back(feature);
  map.getColumnHeaders()[0].filename = "TMT-channel";
  map.getColumnHeaders()[0].label = "126";
  MapPeptideIdentifications ids;
  ids.proteins = {protein()};
  ids.features[{0}] = {peptide()};
  IdentificationDataConverter::importMapIdentifications(map, ids);
  map.eraseMatches([](const auto&, const auto&, const auto& match) { return match.representation != "PEPTIDE"; });
  TEST_REAL_SIMILAR(map[0].getIntensity(), 999)
  TEST_EQUAL(IdentificationDataConverter::exportMapIdentifications(map).at({0})[0].getHits().size(), 1)
  TEST_EQUAL(map.getColumnHeaders()[0].label, "126")
  // identifications of a feature the map does not have are refused, and the map stays as it is
  auto broken = ids;
  broken.features[{500}] = broken.features[{0}];
  TEST_EXCEPTION(Exception::InvalidParameter, IdentificationDataConverter::importMapIdentifications(map, broken))
  TEST_REAL_SIMILAR(map[0].getIntensity(), 999)
  TEST_EQUAL(IdentificationDataConverter::exportMapIdentifications(map).at({0})[0].getHits().size(), 1)
}
END_SECTION

END_TEST
