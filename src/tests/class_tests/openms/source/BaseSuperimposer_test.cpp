// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class TestSuperimposer 
	: public BaseSuperimposer
{
  public:
	TestSuperimposer() 
		: BaseSuperimposer()
	{ 
		check_defaults_ = false; 
	}

	void run(const ConsensusMap& , const ConsensusMap& , TransformationDescription& transformation) override
	{
		Param params;
		params.setValue("slope",1.1);
		params.setValue("intercept", 5.0);
		transformation.fitModel("linear", params);
	}
};

START_TEST(BaseSuperimposer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestSuperimposer* ptr = nullptr;
TestSuperimposer* nullPointer = nullptr;
START_SECTION((BaseSuperimposer()))
	ptr = new TestSuperimposer();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~BaseSuperimposer()))
	delete ptr;
END_SECTION

START_SECTION((virtual void run(const ConsensusMap& map_model, const ConsensusMap& map_scene, TransformationDescription& transformation)=0))
{
  TransformationDescription transformation;
  TestSuperimposer si;
	std::vector<ConsensusMap> maps;
	maps.resize(2);
  si.run(maps[0], maps[1], transformation);
  TEST_STRING_EQUAL(transformation.getModelType(), "linear");
  Param params = transformation.getModelParameters();
  TEST_REAL_SIMILAR(params.getValue("slope"), 1.1)
  TEST_REAL_SIMILAR(params.getValue("intercept"), 5.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



