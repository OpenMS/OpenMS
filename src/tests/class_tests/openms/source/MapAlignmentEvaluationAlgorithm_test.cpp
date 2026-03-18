// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h>
///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmPrecision.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.h>

#include <OpenMS/KERNEL/Feature.h>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
	class MAEA
	 : public MapAlignmentEvaluationAlgorithm
	{
		public:
			void evaluate(const ConsensusMap&, const ConsensusMap&, const double&, const double&, const Peak2D::IntensityType&, const bool use_charge, double& real) override
			{
				bool x = use_charge;
				x=!x;
				real = 1.5;
			}
	};
}

START_TEST(MapAlignmentEvaluation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MAEA* ptr = nullptr;
MAEA* nullPointer = nullptr;
START_SECTION((MapAlignmentEvaluationAlgorithm()))
	ptr = new MAEA();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~MapAlignmentEvaluationAlgorithm()))
	delete ptr;
END_SECTION

START_SECTION((virtual void evaluate(const ConsensusMap &conensus_map_in, const ConsensusMap &consensus_map_gt, const double &rt_dev, const double &mz_dev, const Peak2D::IntensityType &int_dev, const bool use_charge, double &out)=0))
	MAEA maea;
	ConsensusMap map1;
	ConsensusMap map2;
	double rt_dev, mz_dev;
	Peak2D::IntensityType int_dev;
	double real;
	maea.evaluate(map1, map2, rt_dev, mz_dev, int_dev, false, real);
	TEST_EQUAL(real, 1.5)
END_SECTION

START_SECTION((bool isSameHandle(const FeatureHandle &lhs, const FeatureHandle &rhs, const double &rt_dev, const double &mz_dev, const Peak2D::IntensityType &int_dev, const bool use_charge)))
{
	Feature tmp_feature;
	tmp_feature.setRT(100);
	tmp_feature.setMZ(555);
	tmp_feature.setIntensity(200.0f);
	tmp_feature.setCharge(3);
	tmp_feature.setUniqueId(1);

  Feature tmp_feature2;
	tmp_feature2.setRT(101);
	tmp_feature2.setMZ(556);
	tmp_feature2.setIntensity(1199.0f);
	tmp_feature2.setCharge(4);
  tmp_feature2.setUniqueId(2);

	FeatureHandle a(0,tmp_feature);
	FeatureHandle b(0,tmp_feature2);

	MAEA maea;

	TEST_EQUAL(maea.isSameHandle(a, b, 2, 1.5, 1000, false), true);
	TEST_EQUAL(maea.isSameHandle(a, b, 2, 1.5, 1000, true), false);

	tmp_feature2.setCharge(3); // now charge is equal
	FeatureHandle b2(0,tmp_feature2);

	TEST_EQUAL(maea.isSameHandle(a, b2, 2, 1.5, 1000, false), true);
	TEST_EQUAL(maea.isSameHandle(a, b2, 2, 1.5, 1000, true), true);

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
