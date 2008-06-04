// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef FeatureMap<> ElementMapType;

class TestSuperimposer 
	: public BaseSuperimposer<ElementMapType>
{
  public:
	TestSuperimposer() 
		: BaseSuperimposer<ElementMapType>()
	{ 
		check_defaults_ = false; 
	}

	virtual void run(const std::vector<ElementMapType>& maps, TransformationDescription& mapping)
	{
		if (maps.size()!=2) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"scene_map");
		mapping.setName("linear");
		mapping.setParam("slope",1.1);
		mapping.setParam("intercept",5.0);
	}
};

START_TEST(BaseSuperimposer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestSuperimposer* ptr = 0;
CHECK((BaseSuperimposer()))
	ptr = new TestSuperimposer();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~BaseSuperimposer()))
	delete ptr;
RESULT

CHECK(virtual void run(const std::vector<ElementMapType>& maps, TransformationDescription& mapping))
  TransformationDescription mapping;
  TestSuperimposer si;
	std::vector<ElementMapType> maps;
	TEST_EXCEPTION(Exception::IllegalArgument,si.run(maps,mapping))
	maps.resize(2);
  si.run(maps, mapping);
  TEST_STRING_EQUAL(mapping.getName(),"linear");
  TEST_REAL_EQUAL(mapping.getParameters().getValue("slope"),1.1)
  TEST_REAL_EQUAL(mapping.getParameters().getValue("intercept"),5.0)
RESULT

CHECK(void registerChildren())
  NOT_TESTABLE
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



