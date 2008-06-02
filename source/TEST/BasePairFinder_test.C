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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairFinder.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class TestPairFinder 
	: public BasePairFinder
{
  public:
	TestPairFinder() 
		: BasePairFinder()
	{
		check_defaults_ = false; 
	}
};

START_TEST(BasePairFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestPairFinder* ptr = 0;
CHECK((BasePairFinder()))
	ptr = new TestPairFinder();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~BasePairFinder()))
	delete ptr;
RESULT

CHECK((const ConsensusMap& getModelMap() const))
  ConsensusMap map;
  TestPairFinder bpf;
	bpf.setModelMap(0,map);
  TEST_EQUAL(&(bpf.getModelMap()) == &map,true)
RESULT

CHECK((const ConsensusMap& getSceneMap() const))
  ConsensusMap map;
  TestPairFinder bpf;
	bpf.setSceneMap(1,map);
  TEST_EQUAL(&(bpf.getSceneMap()) == &map,true)
RESULT

CHECK((void setModelMap(const ConsensusMap& element_map)))
	NOT_TESTABLE; // see getModelMap()
RESULT

CHECK((void setSceneMap(const ConsensusMap& element_map)))
	NOT_TESTABLE; // see getSceneMap()
RESULT

CHECK(void registerChildren())
  NOT_TESTABLE
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



