// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/METADATA/SampleTreatment.h>
#include <sstream>

///////////////////////////

START_TEST(SampleTreatment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PRECISION(0.001)

class Test: public SampleTreatment
{
	public:
		Test():
		SampleTreatment("Test")
		{
			
		};
		
		Test(const Test& source):
		SampleTreatment(source),
		dummy(source.dummy)
		{
			
		};
		
		virtual ~Test()
		{
			
		};
		
		Test& operator = (const Test& source)
		{
		  if (&source != this)
		  {
		  	SampleTreatment::operator=(source);
		    dummy = source.dummy;
		  }
		  return *this;
		};
	
		
		virtual SampleTreatment* clone() const
		{
			return new Test(*this);
		};

		virtual bool operator== (const SampleTreatment& rhs) const
		{
			return type_==rhs.getType();
		};
	
	public:
		String dummy;
};

// default ctor
Test* dv_ptr = 0;
CHECK([EXTRA] Test())
	dv_ptr = new Test;
	TEST_NOT_EQUAL(dv_ptr, 0)
RESULT

// destructor
CHECK(~SampleTreatment())
	delete dv_ptr;
RESULT

//basic accessors
CHECK([EXTRA] dummy member)
	Test s;
	s.dummy = "TTEST";
	TEST_EQUAL(s.dummy,"TTEST")
RESULT

//getType
CHECK(const String& getType() const)
	Test s;
	TEST_EQUAL(s.getType(),"Test")
RESULT

//meta info
CHECK([EXTRA] MetaInfo)
	Test s;
	//empty
	TEST_EQUAL(s.isMetaEmpty(),true)
	s.metaRegistry().registerName("origin","","");
	s.metaRegistry().registerName("size","","");
	
	s.setMetaValue("origin",string("cow"));
	s.setMetaValue("size",1.0);
	TEST_EQUAL(s.isMetaEmpty(),false)
	TEST_EQUAL(string(s.getMetaValue("origin")),"cow")
	TEST_REAL_EQUAL(double(s.getMetaValue("size")),1.0)
RESULT

//copy ctr
CHECK(SampleTreatment(const SampleTreatment&))
	Test s;
	//set
	s.dummy= "TTEST";
	s.setMetaValue("origin",string("horse"));
	//copy
	Test s2(s);
	//get
	TEST_EQUAL(s2.dummy,"TTEST")
	TEST_EQUAL("horse",s.getMetaValue("origin"))
RESULT

//assignment operator
//copy ctr
CHECK(SampleTreatment(const SampleTreatment&))
	Test s,s2;
	//set
	s.dummy= "TTEST";
	s.setMetaValue("origin",string("horse"));
	//assign
	s2 = s;
	//get
	TEST_EQUAL(s2.dummy,"TTEST")
	TEST_EQUAL("horse",s.getMetaValue("origin"))
RESULT

//clone
CHECK(SampleTreatment& operator=(const SampleTreatment&))
	Test s;
	SampleTreatment* st1;
	SampleTreatment* st;
	Test* dp;
	
	//set
	s.dummy= "TTEST";
	s.setMetaValue("origin",string("horse"));

	//assign
	st1 = &s;
	st = st1->clone();
	dp = dynamic_cast<Test*>(st);
	
	//get
	TEST_EQUAL(dp->dummy,"TTEST")
	TEST_EQUAL("horse",dp->getMetaValue("origin"))
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
