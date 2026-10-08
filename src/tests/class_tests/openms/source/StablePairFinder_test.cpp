// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Feature.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DPosition<2> PositionType;

START_TEST(StablePairFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

StablePairFinder* ptr = nullptr;
StablePairFinder* nullPointer = nullptr;
START_SECTION((StablePairFinder()))
	ptr = new StablePairFinder();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~StablePairFinder()))
	delete ptr;
END_SECTION

START_SECTION((void run(const std::vector<ConsensusMap>& input_maps, ConsensusMap &result_map)))
{
  std::vector<ConsensusMap> input(2);
  Feature feat1;
  Feature feat2;
  Feature feat3;
  PositionType pos1(0,0);
  PositionType pos2(200,300);
  PositionType pos3(400,500);
  feat1.setPosition(pos1);
  feat1.setIntensity(100.0f);
  feat1.setUniqueId(0);
  feat2.setPosition(pos2);
  feat2.setIntensity(300.0f);
  feat2.setUniqueId(1);
  feat3.setPosition(pos3);
  feat3.setIntensity(400.0f);
  feat3.setUniqueId(2);
  ConsensusFeature cons1(0,feat1);
  ConsensusFeature cons2(0,feat2);
  ConsensusFeature cons3(0,feat3);
  input[0].push_back(cons1);
  input[0].push_back(cons2);
  input[0].push_back(cons3);

  Feature feat4;
  Feature feat5;
  Feature feat6;
  PositionType pos4(4,0.04);
  PositionType pos5(204,300.04);
  PositionType pos6(404,500.04);
  feat4.setPosition(pos4);
  feat4.setIntensity(100.0f);
  feat4.setUniqueId(0);
  feat5.setPosition(pos5);
  feat5.setIntensity(300.0f);
  feat5.setUniqueId(1);
  feat6.setPosition(pos6);
  feat6.setIntensity(400.0f);
  feat6.setUniqueId(2);
  ConsensusFeature cons4(1,feat4);
  ConsensusFeature cons5(1,feat5);
  ConsensusFeature cons6(1,feat6);
  input[1].push_back(cons4);
  input[1].push_back(cons5);
  input[1].push_back(cons6);

  input[0].updateRanges();
  input[1].updateRanges();

  StablePairFinder spf;
	Param param = spf.getDefaults();
	spf.setParameters(param);
	ConsensusMap result;
	spf.run(input,result);
	TEST_EQUAL(result.size(),3);
	ABORT_IF(result.size()!=3);

  ConsensusFeature::HandleSetType group1 = result[0].getFeatures();
  ConsensusFeature::HandleSetType group2 = result[1].getFeatures();
  ConsensusFeature::HandleSetType group3 = result[2].getFeatures();

  FeatureHandle ind1(0,feat1);
  FeatureHandle ind2(0,feat2);
  FeatureHandle ind3(0,feat3);
  FeatureHandle ind4(1,feat4);
  FeatureHandle ind5(1,feat5);
  FeatureHandle ind6(1,feat6);

  ConsensusFeature::HandleSetType::const_iterator it;
	it = group1.begin();
  STATUS(*it);
	STATUS(ind1);
	TEST_EQUAL(*(it) == ind1, true)
	++it;
  STATUS(*it);
	STATUS(ind4);
  TEST_EQUAL(*(it) == ind4, true)
	it = group2.begin();
  STATUS(*it);
	STATUS(ind2);
  TEST_EQUAL(*(it) == ind2, true)
	++it;
  STATUS(*it);
	STATUS(ind5);
  TEST_EQUAL(*(it) == ind5, true)
  it = group3.begin();
  STATUS(*it);
	STATUS(ind3);
  TEST_EQUAL(*(it) == ind3, true)
	++it;
  STATUS(*it);
	STATUS(ind6);
  TEST_EQUAL(*(it) == ind6, true)
}
END_SECTION

START_SECTION(([EXTRA] void run(const std::vector<ConsensusMap>& input_maps, ConsensusMap &result_map)))
{
	// test quality calculation:
	std::vector<ConsensusMap> input(2);
  Feature feat1, feat2, feat3;
  PositionType pos1(100, 100);
  PositionType pos2(200, 200);
  PositionType pos3(300, 300);
  feat1.setPosition(pos1);
  feat1.setIntensity(100.0);
  feat1.setUniqueId(0);
  feat2.setPosition(pos2);
  feat2.setIntensity(200.0);
  feat2.setUniqueId(1);
  feat3.setPosition(pos3);
  feat3.setIntensity(300.0);
  feat3.setUniqueId(2);
  // ConsensusFeature cons1(feat1);
  // ConsensusFeature cons2(feat2);
  // ConsensusFeature cons3(feat3);

  StablePairFinder spf;
	Param param = spf.getDefaults();
	param.setValue("distance_RT:max_difference", 1000.0);
	param.setValue("distance_MZ:max_difference", 1000.0);
	param.setValue("second_nearest_gap", 2.0);
	spf.setParameters(param);
	ConsensusMap result;

	// best case:
  input[0].push_back(ConsensusFeature(0, feat1));
	input[1].push_back(ConsensusFeature(1, feat1));

  input[0].updateRanges();
  input[1].updateRanges();

  spf.run(input, result);
	TEST_EQUAL(result.size(), 1);
	TEST_EQUAL(result[0].size(), 2);
	TEST_EQUAL(result[0].getQuality(), 1.0);
	input[0] = result;
	input[1][0] = ConsensusFeature(2, feat1);

  input[0].updateRanges();
  input[1].updateRanges();

	spf.run(input, result);
	TEST_EQUAL(result.size(), 1);
	TEST_EQUAL(result[0].size(), 3);
	TEST_EQUAL(result[0].getQuality(), 1.0);

	// worst case:
	input[0].clear();
	spf.run(input, result);
	TEST_EQUAL(result.size(), 1);
	TEST_EQUAL(result[0].size(), 1);
	TEST_EQUAL(result[0].getQuality(), 0.0);

	// intermediate cases:
	// basis: feat1 < feat2 < feat3
	input[1].clear();
	input[0].push_back(ConsensusFeature(0, feat1));
	input[1].push_back(ConsensusFeature(1, feat2));

  input[0].updateRanges();
  input[1].updateRanges();

	spf.run(input, result);
	ConsensusFeature cons1 = result[0];
	TEST_EQUAL(cons1.size(), 2);
	input[0] = result;
	input[1][0] = ConsensusFeature(2, feat3);

  input[0].updateRanges();
  input[1].updateRanges();

  spf.run(input, result);
	ConsensusFeature cons2 = result[0];
	TEST_EQUAL(cons2.size(), 3);
	TEST_EQUAL(cons1.getQuality() > 0.0, true);
	TEST_EQUAL(cons2.getQuality() > 0.0, true);
	TEST_EQUAL(cons1.getQuality() < 1.0, true);
	TEST_EQUAL(cons2.getQuality() < 1.0, true);
	// quality(feat1, feat2) > quality((feat1, feat2), feat3):
	TEST_EQUAL(cons1.getQuality() > cons2.getQuality(), true);
	input[0].clear();
	input[0].push_back(ConsensusFeature(1, feat2));

  input[0].updateRanges();
  input[1].updateRanges();

  spf.run(input, result);
	ConsensusFeature cons3 = result[0];
	// quality(feat2, feat3) > quality(feat1, feat2), feat3):
	TEST_EQUAL(cons3.getQuality() > cons2.getQuality(), true);
	input[0][0] = ConsensusFeature(0, feat1);

  input[0].updateRanges();
  input[1].updateRanges();

  spf.run(input, result);
	ConsensusFeature cons4 = result[0];
	// quality(feat1, feat3) < quality(feat1, feat2), feat3):
	TEST_EQUAL(cons4.getQuality() < cons2.getQuality(), true);
}
END_SECTION

START_SECTION(([EXTRA] run() with use_identifications on identification data))
{
  using ID = IdentificationData;
  std::vector<ConsensusMap> input(2);
  const auto annotate = [](ConsensusMap& map, ConsensusFeature& feature, const std::string& sequence) {
    auto& data = map.getIdentificationData();
    if (data.getRuns().empty())
    {
      auto& run = data.addRun("search");
      ID::ScoreDefinition score;
      score.name = "score";
      score.higher_better = true;
      run.setPrimaryScore(run.addScore(score));
      run.addSource({});
    }
    auto& run = data.getRun("search");
    ID::MatchData match;
    match.representation = sequence;
    const auto query = run.addIdentification(run.getSources()[0].id, ID::Observation {});
    feature.addIDMatch({run.getUuid(), run.addMatch(query, match, {1.0})});
  };
  Feature a, b, c;
  a.setPosition(PositionType(0, 0));
  b.setPosition(PositionType(4, 0.04));
  c.setPosition(PositionType(1, 0.01)); // nearer, but with another peptide
  for (Feature* feature : {&a, &b, &c})
  {
    feature->setIntensity(100.0f);
  }
  a.setUniqueId(0);
  b.setUniqueId(0);
  c.setUniqueId(1);
  ConsensusFeature cons_a(0, a), cons_b(1, b), cons_c(1, c);
  annotate(input[0], cons_a, "PEPA");
  annotate(input[1], cons_b, "PEPA");
  annotate(input[1], cons_c, "PEPB");
  input[0].push_back(cons_a);
  input[1].push_back(cons_b);
  input[1].push_back(cons_c);
  input[0].updateRanges();
  input[1].updateRanges();

  // the partner of feature "a" in the second map:
  const auto partner = [](const ConsensusMap& result) -> Int {
    for (const auto& feature : result)
    {
      if (feature.size() == 2) return Int(feature.rbegin()->getUniqueId());
    }
    return -1;
  };
  StablePairFinder spf;
  Param param = spf.getDefaults();
  spf.setParameters(param);
  ConsensusMap result;
  spf.run(input, result);
  TEST_EQUAL(result.size(), 2)
  TEST_EQUAL(partner(result), 1)

  param.setValue("use_identifications", "true");
  spf.setParameters(param);
  spf.run(input, result);
  TEST_EQUAL(result.size(), 2)
  TEST_EQUAL(partner(result), 0)

  // the result has the identifications of both maps:
  const auto& data = result.getIdentificationData();
  TEST_EQUAL(data.getRuns().size(), 2)
  for (const auto& feature : result)
  {
    TEST_EQUAL(feature.getLinkedIdentifications(data).size(), feature.size())
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
