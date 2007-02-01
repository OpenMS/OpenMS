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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>

#include <OpenMS/CONCEPT/Exception.h>

///////////////////////////

START_TEST(BaseSeeder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

class TestSeeder : public BaseSeeder
{
  public:
	TestSeeder()
		: BaseSeeder()
	{
		setName(TestSeeder::getName());
		check_defaults_ = false;
	}
	
	IndexSet nextSeed() throw (NoSuccessor)
	{
		FeaFiModule::IndexSet  set;
		set.insert(std::make_pair(7,7));
		return set;
	}
	
	static const String getName(){ return "TestSeeder"; }
};

// default ctor
TestSeeder* ptr = 0;
CHECK(TestSeeder())
	ptr = new TestSeeder();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(~TestSeeder())
	delete ptr;
RESULT

// assignment operator
CHECK(TestSeeder& operator = (const TestSeeder& source))
	TestSeeder tm1;
  TestSeeder tm2;
  tm2 = tm1;

  TestSeeder tm3;

  tm1 = TestSeeder();
	TEST_EQUAL(tm3,tm2)
RESULT

// copy constructor
CHECK(TestSeeder(const TestSeeder& source))
	TestSeeder fp1;	
 
  TestSeeder fp2(fp1);

  TestSeeder fp3;
 
  fp1 = TestSeeder();
	TEST_EQUAL(fp2, fp3)
RESULT

CHECK(IndexSet nextSeed() throw (NoSuccessor))
	TestSeeder s;
	FeaFiModule::IndexSet  almost_empty;
	almost_empty.insert(std::make_pair(7,7));
	TEST_EQUAL(s.nextSeed()==almost_empty,true);
RESULT

CHECK(static void registerChildren())
	// not much happening here
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
