// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

namespace
{
  using ID = IdentificationData;

  /// A map with a feature at (@p rt, 500) that links a match of "PEPTIDE", and an unassigned identification
  FeatureMap makeMap(double rt, UInt64 unique_id)
  {
    FeatureMap map;
    auto& run = map.getIdentificationData().addRun("search");
    ID::ScoreDefinition score;
    score.name = "score";
    score.higher_better = true;
    run.setPrimaryScore(run.addScore(score));
    const auto source = run.addSource({});
    ID::MatchData match;
    match.representation = "PEPTIDE";
    Feature feature;
    feature.setRT(rt);
    feature.setMZ(500.0);
    feature.setIntensity(100.0f);
    feature.setUniqueId(unique_id);
    feature.addIDMatch({run.getUuid(), run.addMatch(run.addIdentification(source, ID::Observation {}), match, {1.0})});
    run.addIdentification(source, ID::Observation {});
    map.push_back(feature);
    map.updateRanges();
    return map;
  }

  /// Whether the identifications of each run of @p data have the run's position as map index
  bool markedByRun(const ID& data)
  {
    for (Size i = 0; i < data.getRuns().size(); ++i)
    {
      for (const auto& query : data.getRuns()[i].getSources()[0].identifications)
      {
        if (Size(query.getMetaValue("map_index", -1)) != i) return false;
      }
    }
    return true;
  }
} // namespace


START_TEST(FeatureGroupingAlgorithmUnlabeled, "$Id FeatureGroupingAlgorithmUnlabeled_test.C 139 2006-07-14 10:08:39Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureGroupingAlgorithmUnlabeled* ptr = nullptr;
FeatureGroupingAlgorithmUnlabeled* nullPointer = nullptr;
START_SECTION((FeatureGroupingAlgorithmUnlabeled()))
	ptr = new FeatureGroupingAlgorithmUnlabeled();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithmUnlabeled()))
	delete ptr;
END_SECTION

START_SECTION((virtual void group(const std::vector< FeatureMap > &maps, ConsensusMap &out)))
	// This is tested extensively in TEST/TOPP
	NOT_TESTABLE;
END_SECTION

START_SECTION(([EXTRA] group() with identification data))
{
  // the third map has the runs of the second (they stay apart):
  std::vector<FeatureMap> maps {makeMap(100.0, 1), makeMap(101.0, 2)};
  maps.push_back(maps[1]);
  maps[2][0].setRT(102.0);
  maps[2].updateRanges();
  FeatureGroupingAlgorithmUnlabeled algo;
  ConsensusMap out;
  algo.group(maps, out);
  TEST_EQUAL(out.size(), 1)
  ABORT_IF(out.size() != 1)
  TEST_EQUAL(out[0].size(), 3)
  const auto& data = out.getIdentificationData();
  TEST_EQUAL(data.getRuns().size(), 3)
  TEST_EQUAL(markedByRun(data), true)
  // each run keeps its unassigned identification:
  for (const auto& run : data.getRuns())
  {
    TEST_EQUAL(run.getNumberOfIdentifications(), 2)
  }
  TEST_EQUAL(out[0].getLinkedIdentifications(data).size(), 3)
  TEST_EQUAL(out.getIdentificationData().getRuns()[0].getUuid(), maps[0].getIdentificationData().getRuns()[0].getUuid())
  TEST_NOT_EQUAL(data.getRuns()[2].getUuid(), data.getRuns()[1].getUuid())
}
END_SECTION

START_SECTION((void addToGroup(int map_id, const FeatureMap& feature_map)))
{
  std::vector<FeatureMap> maps {makeMap(100.0, 1), makeMap(101.0, 2)};
  FeatureGroupingAlgorithmUnlabeled algo;
  algo.setReference(0, maps[0]);
  algo.addToGroup(1, maps[1]);
  ConsensusMap result = algo.getResultMap();
  TEST_EQUAL(result.size(), 1)
  ABORT_IF(result.size() != 1)
  TEST_EQUAL(result[0].size(), 2)
  // the raw result links the identifications of both maps:
  TEST_EQUAL(result.getIdentificationData().getRuns().size(), 2)
  TEST_EQUAL(result[0].getLinkedIdentifications(result.getIdentificationData()).size(), 2)
  // they get the map index, like after group():
  FeatureGroupingAlgorithm::groupIdentifications(maps, result);
  TEST_EQUAL(markedByRun(result.getIdentificationData()), true)
  TEST_EQUAL(result[0].getLinkedIdentifications(result.getIdentificationData()).size(), 2)
}
END_SECTION

START_SECTION(([EXTRA] group() with peptide identifications))
{
  // maps with peptide identifications are grouped on identification data; so does the result have them:
  std::vector<FeatureMap> maps {makeMap(100.0, 1), makeMap(101.0, 2)};
  for (auto& map : maps)
  {
    IdentificationDataConverter::exportFeatureIDs(map);
  }
  TEST_EQUAL(maps[0][0].getPeptideIdentifications().size(), 1)
  FeatureGroupingAlgorithmUnlabeled algo;
  ConsensusMap out;
  algo.group(maps, out);
  TEST_EQUAL(out.size(), 1)
  ABORT_IF(out.size() != 1)
  TEST_EQUAL(out.getIdentificationData().empty(), true)
  TEST_EQUAL(out.getProteinIdentifications().size(), 2)
  TEST_EQUAL(out[0].getPeptideIdentifications().size(), 2)
  ABORT_IF(out[0].getPeptideIdentifications().size() != 2)
  TEST_EQUAL(Size(out[0].getPeptideIdentifications()[0].getMetaValue("map_index")), 0)
  TEST_EQUAL(Size(out[0].getPeptideIdentifications()[1].getMetaValue("map_index")), 1)
  TEST_EQUAL(out.getUnassignedPeptideIdentifications().size(), 2)
  ABORT_IF(out.getUnassignedPeptideIdentifications().size() != 2)
  TEST_EQUAL(Size(out.getUnassignedPeptideIdentifications()[1].getMetaValue("map_index")), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



