// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

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

	virtual void run(const ConsensusMap& map_model, const ConsensusMap& map_scene, TransformationDescription& transformation)
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

TestSuperimposer* ptr = 0;
TestSuperimposer* nullPointer = 0;
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
  Param params;
	transformation.getModelParameters(params);
  TEST_REAL_SIMILAR(params.getValue("slope"), 1.1)
  TEST_REAL_SIMILAR(params.getValue("intercept"), 5.0)
}
END_SECTION

START_SECTION((static void registerChildren()))
  NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



