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
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class TestPairFinder 
	: public BaseGroupFinder
{
  public:
	TestPairFinder() 
		: BaseGroupFinder()
	{
		check_defaults_ = false; 
	}
	void run(const std::vector<ConsensusMap>&, ConsensusMap&) override
	{
	}
};

START_TEST(BaseGroupFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestPairFinder* ptr = nullptr;
TestPairFinder* nullPointer = nullptr;
START_SECTION((BaseGroupFinder()))
	ptr = new TestPairFinder();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~BaseGroupFinder()))
	delete ptr;
END_SECTION


START_SECTION((virtual void run(const std::vector< ConsensusMap > &input, ConsensusMap &result)=0))
	NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



