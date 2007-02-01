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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/CONCEPT/Exception.h>


///////////////////////////

START_TEST(BaseExtender, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

class TestExtender : public BaseExtender
{
  public:
	TestExtender(): BaseExtender()
	{
		setName(TestExtender::getName());
		check_defaults_ = false;
	}

	const IndexSet& extend(const IndexSet& /*seed_region*/)
	{
		return region_;	
	}
	
	static const String getName(){ return "TestExtender"; }

};

// default ctor
TestExtender* ptr = 0;
CHECK(TestExtender())
	ptr = new TestExtender();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(~TestExtender())
	delete ptr;
RESULT

// assignment operator
CHECK(TestExtender& operator = (const TestExtender& source))
	TestExtender tm1;
  TestExtender tm2;
  tm2 = tm1;

  TestExtender tm3;

  tm1 = TestExtender();
	TEST_EQUAL(tm3,tm2)
RESULT

// copy constructor
CHECK(TestExtender(const TestExtender& source))
	TestExtender fp1;	
 
  TestExtender fp2(fp1);

  TestExtender fp3;
 
  fp1 = TestExtender();
	TEST_EQUAL(fp2, fp3)
RESULT

CHECK(const IndexSet& extend(const IndexSet& seed_region))
	TestExtender text;
	FeaFiModule::IndexSet  inds;
	inds.insert(std::make_pair(7,7));
	FeaFiModule::IndexSet result = text.extend(inds);
  TEST_EQUAL(result.size(),0)
RESULT

CHECK(static void registerChildren())
	// not much happening here
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
