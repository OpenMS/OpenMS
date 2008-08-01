// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

#include <vector>

///////////////////////////

START_TEST(TransformationDescription, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


TransformationDescription* ptr = 0;
CHECK((TransformationDescription()))
	ptr = new TransformationDescription;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~TransformationDescription()))
	delete ptr;
RESULT

CHECK((const String& getName() const))
	TransformationDescription td;
	TEST_STRING_EQUAL(td.getName(),"")
RESULT

CHECK((void setName(const String& name)))
	TransformationDescription td;
	td.setName("bla");
	TEST_STRING_EQUAL(td.getName(),"bla")

RESULT

CHECK((const Param& getParameters() const))
	TransformationDescription td;
	TEST_EQUAL(td.getParameters(),Param())
RESULT

CHECK((DoubleReal getParam(const String& name) const))
	TransformationDescription td;
	TEST_EXCEPTION(Exception::ElementNotFound, td.getParam("bla"))
RESULT

CHECK((void setParam(const String& name, DoubleReal value)))
	TransformationDescription td;
	td.setParam("bla",4.5);
	TEST_REAL_EQUAL(td.getParam("bla"),4.5)
RESULT

CHECK((void setParameters(const Param& param)))
	TransformationDescription td;
	Param p;
	p.setValue("int",5);
	td.setParameters(p);
	TEST_EQUAL((Int)td.getParameters().size(),1)
	TEST_EQUAL((Int)td.getParameters().getValue("int"),5)
RESULT


TransformationDescription::PairVector pairs;
pairs.push_back(make_pair(1.2,5.2));
pairs.push_back(make_pair(2.2,6.25));
pairs.push_back(make_pair(3.2,7.3));

CHECK(const PairVector& getPairs() const)
	TransformationDescription td;
	TEST_EQUAL(td.getPairs().size(),0)
RESULT

CHECK(void setPairs(const PairVector& pairs))	
{
	TransformationDescription td;
	td.setPairs(pairs);
	TEST_EQUAL(td.getPairs().size(),3)

	TransformationDescription::PairVector pairs_empty;
  pairs_empty.clear();
	td.setPairs(pairs_empty);
	TEST_EQUAL(td.getPairs().size(),0)
}
RESULT

CHECK((TransformationDescription(const TransformationDescription& rhs)))
{
	TransformationDescription td;
	td.setName("dummy");
	td.setParam("int",5);
	td.setPairs(pairs);
	
	TEST_EQUAL(td.getName()==td.getName(),true)
	TEST_EQUAL(td.getParameters()==td.getParameters(),true)
	TEST_EQUAL(td.getPairs().size(),3)
}
RESULT

CHECK((TransformationDescription& operator = (const TransformationDescription& rhs)))
{
	TransformationDescription td;
	td.setName("dummy");
	td.setParam("int",5);
	td.setPairs(pairs);
	TransformationDescription td2;
	td2 = td;
	
 	TEST_STRING_EQUAL(td2.getName(),td.getName());
	TEST_EQUAL(td2.getParameters()==td.getParameters(),true);
	TEST_EQUAL(td2.getPairs()==td.getPairs(),true);
}
RESULT

CHECK((void clear()))
{
	TransformationDescription td;

	td.setName("linear");
	td.setParam("slope",2.0);
	td.setParam("intercept",47.12);
	td.setPairs(pairs);

	DoubleReal value = 5.0;
	td.apply(value);
	TEST_REAL_EQUAL(value,57.12);

 	TEST_STRING_EQUAL(td.getName(),"linear");
	TEST_EQUAL((DoubleReal)td.getParameters().getValue("slope"),2.0);
	TEST_EQUAL((DoubleReal)td.getParameters().getValue("intercept"),47.12);
	TEST_EQUAL(td.getPairs()==pairs,true);
	TEST_EQUAL(td.getPairs().size(),3);

	td.clear();

 	TEST_STRING_EQUAL(td.getName(),"");
	TEST_EQUAL(td.getParameters().empty(),true);
	TEST_EQUAL(td.getPairs()==pairs,false);
	TEST_EQUAL(td.getPairs().size(),0);
	TEST_EXCEPTION(Exception::IllegalArgument,td.apply(value));
}
RESULT

CHECK((void apply(DoubleReal& value)))
	DoubleReal value = 5.0;
	TransformationDescription td;
	
	//test missing name and parameters
 	TEST_EXCEPTION(Exception::IllegalArgument,td.apply(value))

	td.setName("bla");
	TEST_EXCEPTION(Exception::IllegalArgument,td.apply(value))
	
	//test with identity
	td.setName("none");
	td.apply(value);
	TEST_REAL_EQUAL(value,5.0);
	
	//test for missing parameter
	td.setName("linear");
	td.setParam("slope",1.0);	
	TEST_EXCEPTION(Exception::IllegalArgument,td.apply(value))
	
	//real test (linear, identity)
	td.setParam("intercept",0.0);
	TEST_REAL_EQUAL(value,5.0);

	//real test (linear, no identity)
	td.setParam("slope",2.0);
	td.setParam("intercept",47.12);
	td.apply(value);
	TEST_REAL_EQUAL(value,57.12);
		
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
