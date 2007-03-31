// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModelFitter.h>
#include <OpenMS/CONCEPT/Exception.h>


///////////////////////////

START_TEST(BaseModelFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

class TestFitter : public BaseModelFitter
{
  public:
	TestFitter(): BaseModelFitter()
	{
		setName(TestFitter::getName());
		check_defaults_ = false;
	}

	Feature fit(const ChargedIndexSet& /*extension*/) throw(UnableToFit)
	{
		Feature f;
		return f;	
	}
	
	static const String getName(){ return "TestFitter"; }
};

// default ctor
TestFitter* ptr = 0;
CHECK((BaseModelFitter()))
	ptr = new TestFitter();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK((virtual ~BaseModelFitter()))
	delete ptr;
RESULT

// assignment operator
CHECK((virtual BaseModelFitter& operator=(const BaseModelFitter &source)))
	TestFitter tm1;
  TestFitter tm2;
  tm2 = tm1;

  TestFitter tm3;

  tm1 = TestFitter();
	TEST_EQUAL(tm3,tm2)
RESULT

// copy constructor
CHECK((BaseModelFitter(const BaseModelFitter &source)))
	TestFitter fp1;	
 
  TestFitter fp2(fp1);

  TestFitter fp3;
 
  fp1 = TestFitter();
	TEST_EQUAL(fp2, fp3)
RESULT

CHECK((virtual Feature fit(const ChargedIndexSet &extension)=0 throw (UnableToFit)))
	TestFitter ft;
	FeaFiModule::ChargedIndexSet  inds;
	inds.insert(std::make_pair(7,7));
	Feature result = ft.fit(inds);
	Feature empty;
  TEST_EQUAL(result==empty,true)
RESULT
	
CHECK((static void registerChildren()))
	// not much happening here
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
