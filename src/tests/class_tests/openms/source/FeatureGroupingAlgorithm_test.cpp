// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>
///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
	class FGA
	 : public FeatureGroupingAlgorithm
	{
		public:
			void group(const vector< FeatureMap >&, ConsensusMap& map) override
			{
			  map.getColumnHeaders()[0].filename = "bla";
				map.getColumnHeaders()[0].size = 5;
			}
	};
}

namespace
{
  using ID = IdentificationData;

  /// A peptide search run (score: higher is better) with one source
  ID::Run& addSearchRun(ID& data, const std::string& name)
  {
    auto& run = data.addRun(name);
    ID::ScoreDefinition score;
    score.name = "score";
    score.higher_better = true;
    run.setPrimaryScore(run.addScore(score));
    run.addSource({});
    return run;
  }

  ID::QueryId addQuery(ID::Run& run)
  {
    return run.addIdentification(run.getSources()[0].id, ID::Observation {});
  }

  ID::MatchId addPeptide(ID::Run& run, ID::QueryId query, const std::string& sequence)
  {
    ID::MatchData match;
    match.representation = sequence;
    return run.addMatch(query, match, {1.0});
  }
} // namespace

START_TEST(FeatureGroupingAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FGA* ptr = nullptr;
FGA* nullPointer = nullptr;
START_SECTION((FeatureGroupingAlgorithm()))
	ptr = new FGA();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithm()))
	delete ptr;
END_SECTION

START_SECTION((virtual void group(const vector< FeatureMap > &maps, ConsensusMap &out)=0))
	FGA fga;
	vector< FeatureMap > in;
	ConsensusMap map;
	fga.group(in,map);
	TEST_EQUAL(map.getColumnHeaders()[0].filename, "bla")
END_SECTION

START_SECTION((void transferSubelements(const vector<ConsensusMap>& maps, ConsensusMap& out) const))
{
	vector<ConsensusMap> maps(2);
	maps[0].getColumnHeaders()[0].filename = "file1";
	maps[0].getColumnHeaders()[0].size = 1;
	maps[0].getColumnHeaders()[1].filename = "file2";
	maps[0].getColumnHeaders()[1].size = 1;
	maps[1].getColumnHeaders()[0].filename = "file3";
	maps[1].getColumnHeaders()[0].size = 1;
	maps[1].getColumnHeaders()[1].filename = "file4";
	maps[1].getColumnHeaders()[1].size = 1;

	Feature feat1, feat2, feat3, feat4;

  FeatureHandle handle1(0, feat1), handle2(1, feat2), handle3(0, feat3),
		handle4(1, feat4);

	maps[0].resize(1);
	maps[0][0].insert(handle1);
	maps[0][0].insert(handle2);
	maps[0][0].setUniqueId(1);
	maps[1].resize(1);
	maps[1][0].insert(handle3);
	maps[1][0].insert(handle4);
	maps[1][0].setUniqueId(2);

	ConsensusMap out;
	FeatureHandle handle5(0, static_cast<BaseFeature>(maps[0][0]));
	FeatureHandle handle6(1, static_cast<BaseFeature>(maps[1][0]));
	out.resize(1);
	out[0].insert(handle5);
	out[0].insert(handle6);

	// need an instance of FeatureGroupingAlgorithm:
	FeatureGroupingAlgorithm* algo = new FeatureGroupingAlgorithmKD();

	// identifications (in identification data) whose old map index was saved, or not:
	auto& run = addSearchRun(out.getIdentificationData(), "search");
	const auto saved = addQuery(run), unsaved = addQuery(run), none = addQuery(run);
	ID::Observation observation;
	observation.setMetaValue("map_index", 1); // the second input map ...
	observation.setMetaValue("old_map_index", 1); // ... and its second map
	run.replaceObservation(saved, observation);
	observation.removeMetaValue("old_map_index");
	run.replaceObservation(unsaved, observation);

	algo->transferSubelements(maps, out);

	TEST_EQUAL(out.getColumnHeaders().size(), 4);
	TEST_EQUAL(out.getColumnHeaders()[0].filename, "file1");
	TEST_EQUAL(out.getColumnHeaders()[3].filename, "file4");
	TEST_EQUAL(out.size(), 1);
	TEST_EQUAL(out[0].size(), 4);

	ConsensusFeature::HandleSetType group = out[0].getFeatures();
	ConsensusFeature::HandleSetType::const_iterator it = group.begin();
	handle3.setMapIndex(2);
	handle4.setMapIndex(3);
	TEST_EQUAL(*it++ == handle1, true);
	TEST_EQUAL(*it++ == handle2, true);
	TEST_EQUAL(*it++ == handle3, true);
	TEST_EQUAL(*it++ == handle4, true);

	// identification data: the map index of an identification follows its old one, if it was saved
	TEST_EQUAL(Size(run.getIdentification(saved).getMetaValue("map_index")), 3)
	TEST_EQUAL(run.getIdentification(saved).metaValueExists("old_map_index"), false)
	TEST_EQUAL(run.getIdentification(unsaved).metaValueExists("map_index"), false)
	TEST_EQUAL(run.getIdentification(none).metaValueExists("map_index"), false)
	delete algo;
}
END_SECTION

START_SECTION((static MapIdentifications getMapIdentifications(const FeatureMap& map)))
{
  FeatureMap map;
  auto& run = addSearchRun(map.getIdentificationData(), "search");
  const auto first = addQuery(run), second = addQuery(run), third = addQuery(run);
  const auto match = addPeptide(run, first, "PEPTIDE");
  const auto sub_match = addPeptide(run, second, "PEPTIDER");
  addPeptide(run, third, "PEPTIDEK");
  Feature feature, subordinate;
  feature.addIDMatch({run.getUuid(), match});
  subordinate.addIDMatch({run.getUuid(), sub_match});
  subordinate.addIDQuery({run.getUuid(), third});
  feature.getSubordinates().push_back(subordinate);
  map.push_back(feature);

  const auto identifications = FeatureGroupingAlgorithm::getMapIdentifications(map);
  TEST_EQUAL(identifications.data == map.getIdentificationData(), true)
  TEST_EQUAL(identifications.matches.size(), 2)
  TEST_EQUAL(identifications.matches.contains({run.getUuid(), match}), true)
  TEST_EQUAL(identifications.matches.contains({run.getUuid(), sub_match}), true)
  TEST_EQUAL(identifications.queries.size(), 1)
  TEST_EQUAL(identifications.queries.contains({run.getUuid(), third}), true)
}
END_SECTION

START_SECTION((static MapIdentifications getMapIdentifications(const ConsensusMap& map)))
{
  ConsensusMap map;
  auto& run = addSearchRun(map.getIdentificationData(), "search");
  const auto query = addQuery(run);
  ConsensusFeature feature;
  feature.addIDQuery({run.getUuid(), query});
  map.push_back(feature);
  const auto identifications = FeatureGroupingAlgorithm::getMapIdentifications(map);
  TEST_EQUAL(identifications.queries.size(), 1)
  TEST_EQUAL(identifications.matches.empty(), true)
}
END_SECTION

START_SECTION((static void groupIdentifications(const std::vector<FeatureMap>& maps, ConsensusMap& grouped)))
{
  // first map: a grouped feature, its subordinate and a feature left out
  FeatureMap map0;
  auto& run0 = addSearchRun(map0.getIdentificationData(), "first");
  const auto grouped_query = addQuery(run0), sub_query = addQuery(run0), partial_query = addQuery(run0),
             unassigned_query = addQuery(run0), left_out_query = addQuery(run0);
  const auto grouped_match = addPeptide(run0, grouped_query, "PEPTIDE");
  const auto sub_match = addPeptide(run0, sub_query, "PEPTIDER");
  const auto left_out_match = addPeptide(run0, partial_query, "PEPTIDEK");
  const auto unassigned_match = addPeptide(run0, partial_query, "PEPTIDEH");
  Feature grouped, subordinate, left_out;
  grouped.setUniqueId(1);
  grouped.addIDMatch({run0.getUuid(), grouped_match});
  subordinate.addIDMatch({run0.getUuid(), sub_match});
  grouped.getSubordinates().push_back(subordinate);
  left_out.setUniqueId(2);
  left_out.addIDMatch({run0.getUuid(), left_out_match});
  left_out.addIDQuery({run0.getUuid(), left_out_query});
  map0.push_back(grouped);
  map0.push_back(left_out);
  // second map: a grouped feature
  FeatureMap map1;
  auto& run1 = addSearchRun(map1.getIdentificationData(), "second");
  const auto second_query = addQuery(run1);
  Feature second;
  second.setUniqueId(3);
  second.addIDQuery({run1.getUuid(), second_query});
  map1.push_back(second);

  ConsensusMap out;
  ConsensusFeature consensus;
  consensus.insert(0, map0[0]);
  consensus.insert(1, map1[0]);
  out.push_back(consensus);
  FeatureGroupingAlgorithm::groupIdentifications(std::vector<FeatureMap> {map0, map1}, out);

  const auto& data = out.getIdentificationData();
  TEST_EQUAL(data.getRuns().size(), 2)
  ABORT_IF(data.getRuns().size() != 2)
  const auto& first = data.getRuns()[0];
  TEST_EQUAL(first.getUuid(), run0.getUuid())
  // what only the subordinate or the feature left out links is dropped, unassigned identifications and matches stay:
  TEST_EQUAL(first.getNumberOfIdentifications(), 3)
  TEST_EQUAL(first.findIdentification(grouped_query) != nullptr, true)
  TEST_EQUAL(first.findIdentification(sub_query) == nullptr, true)
  TEST_EQUAL(first.findIdentification(left_out_query) == nullptr, true)
  TEST_EQUAL(first.findIdentification(unassigned_query) != nullptr, true)
  TEST_EQUAL(first.findMatch(left_out_match) == nullptr, true)
  TEST_EQUAL(first.findMatch(unassigned_match) != nullptr, true)
  // the identifications are marked with the index of their map:
  for (const auto& query : first.getSources()[0].identifications)
  {
    TEST_EQUAL(Size(query.getMetaValue("map_index")), 0)
  }
  TEST_EQUAL(Size(data.getRuns()[1].getIdentification(second_query).getMetaValue("map_index")), 1)
  // the grouped feature links both identifications:
  const auto linked = out[0].getLinkedIdentifications(data);
  TEST_EQUAL(linked.size(), 2)
  ABORT_IF(linked.size() != 2)
  TEST_EQUAL(linked[0].query->getId() == grouped_query, true)
  TEST_EQUAL(linked[0].matches.size(), 1)
  TEST_EQUAL(linked[1].query->getId() == second_query, true)
}
END_SECTION

START_SECTION((static void groupIdentifications(std::vector<MapIdentifications> maps, ConsensusMap& grouped)))
{
  // the identifications of a map that its features do not link are kept (and marked):
  FeatureGroupingAlgorithm::MapIdentifications identifications;
  auto& run = addSearchRun(identifications.data, "search");
  const auto query = addQuery(run);
  ConsensusMap out;
  addSearchRun(out.getIdentificationData(), "replaced");
  // (the second map; the first has none)
  std::vector<FeatureGroupingAlgorithm::MapIdentifications> maps(1);
  maps.push_back(identifications);
  FeatureGroupingAlgorithm::groupIdentifications(std::move(maps), out);
  const auto& data = out.getIdentificationData();
  TEST_EQUAL(data.getRuns().size(), 1)
  ABORT_IF(data.getRuns().size() != 1)
  TEST_EQUAL(Size(data.getRuns()[0].getIdentification(query).getMetaValue("map_index")), 1)
}
END_SECTION

START_SECTION((static void groupIdentifications(const std::vector<ConsensusMap>& maps, ConsensusMap& grouped)))
{
  // the identifications of consensus maps get the index of their map in the grouping:
  ConsensusMap map;
  auto& run = addSearchRun(map.getIdentificationData(), "search");
  const auto query = addQuery(run);
  ID::Observation observation;
  observation.setMetaValue("map_index", 5);
  run.replaceObservation(query, observation);
  ConsensusFeature feature;
  feature.setUniqueId(1);
  feature.addIDQuery({run.getUuid(), query});
  map.push_back(feature);
  ConsensusMap out;
  ConsensusFeature grouped;
  grouped.insert(0, map[0]);
  out.push_back(grouped);
  FeatureGroupingAlgorithm::groupIdentifications(std::vector<ConsensusMap> {map}, out);
  TEST_EQUAL(Size(out.getIdentificationData().getRuns()[0].getIdentification(query).getMetaValue("map_index")), 0)
}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
