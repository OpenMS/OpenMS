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
	delete algo;
}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
