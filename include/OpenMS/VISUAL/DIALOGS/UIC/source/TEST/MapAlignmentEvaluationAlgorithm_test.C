// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h>
///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmPrecision.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.h>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
	class MAEA
	 : public MapAlignmentEvaluationAlgorithm
	{
		public:
			void evaluate(const ConsensusMap&, const ConsensusMap&, const DoubleReal&, const DoubleReal&, const Peak2D::IntensityType&, const bool use_charge, DoubleReal& real)
			{
				bool x = use_charge;
				x=!x;
				real = 1.5;
			}
	};
}

START_TEST(MapAlignmentEvaluation, "$Id MapAlignmentEvaluationAlgorithm_test.C 139 2006-07-14 10:08:39Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MAEA* ptr = 0;
START_SECTION((MapAlignmentEvaluationAlgorithm()))
	ptr = new MAEA();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~MapAlignmentEvaluationAlgorithm()))
	delete ptr;
END_SECTION

START_SECTION((virtual void evaluate(const ConsensusMap &conensus_map_in, const ConsensusMap &consensus_map_gt, const DoubleReal &rt_dev, const DoubleReal &mz_dev, const Peak2D::IntensityType &int_dev, const bool use_charge, DoubleReal &out)=0))
	MAEA maea;
	ConsensusMap map1;
	ConsensusMap map2;
	DoubleReal rt_dev, mz_dev;
	Peak2D::IntensityType int_dev;
	DoubleReal real;
	maea.evaluate(map1, map2, rt_dev, mz_dev, int_dev, false, real);
	TEST_EQUAL(real, 1.5)
END_SECTION

START_SECTION((bool isSameHandle(const FeatureHandle &lhs, const FeatureHandle &rhs, const DoubleReal &rt_dev, const DoubleReal &mz_dev, const Peak2D::IntensityType &int_dev, const bool use_charge)))
{
	Feature tmp_feature;
	tmp_feature.setRT(100);
	tmp_feature.setMZ(555);
	tmp_feature.setIntensity(200.0f);
	tmp_feature.setCharge(3);

  Feature tmp_feature2;
	tmp_feature2.setRT(101);
	tmp_feature2.setMZ(556);
	tmp_feature2.setIntensity(1199.0f);
	tmp_feature2.setCharge(4);

	FeatureHandle a(0,1, tmp_feature);
	FeatureHandle b(0,2, tmp_feature2);

	MAEA maea;

	TEST_EQUAL(maea.isSameHandle(a, b, 2, 1.5, 1000, false), true);
	TEST_EQUAL(maea.isSameHandle(a, b, 2, 1.5, 1000, true), false);

	tmp_feature2.setCharge(3); // now charge is equal
	FeatureHandle b2(0,1, tmp_feature2);

	TEST_EQUAL(maea.isSameHandle(a, b2, 2, 1.5, 1000, false), true);
	TEST_EQUAL(maea.isSameHandle(a, b2, 2, 1.5, 1000, true), true);

}
END_SECTION

START_SECTION((static void registerChildren()))
{
	TEST_STRING_EQUAL(Factory<MapAlignmentEvaluationAlgorithm>::registeredProducts()[0],MapAlignmentEvaluationAlgorithmPrecision::getProductName());
	TEST_STRING_EQUAL(Factory<MapAlignmentEvaluationAlgorithm>::registeredProducts()[1],MapAlignmentEvaluationAlgorithmRecall::getProductName());
	TEST_EQUAL(Factory<MapAlignmentEvaluationAlgorithm>::registeredProducts().size(),2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
