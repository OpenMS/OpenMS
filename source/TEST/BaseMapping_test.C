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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseMapping.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class TestMapping : public BaseMapping
{
  public:
  TestMapping() : BaseMapping(){}
  TestMapping(const TestMapping& tm) : BaseMapping(tm){}
  TestMapping& operator=(const TestMapping& tm)
  {
     BaseMapping::operator=(tm);
     return *this;
  }
    
   virtual void apply(DPosition<1>& ) const {}

   virtual void apply(DoubleReal& ) const {};
  
   virtual const String getName() { return "";}
};

START_TEST(BaseMapping, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestMapping* ptr = 0;
CHECK((BaseMapping()))
	ptr = new TestMapping();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~BaseMapping()))
	delete ptr;
RESULT

CHECK((BaseMapping& operator = (const BaseMapping& rhs)))
  TestMapping tm;
  Param param;
  param.setValue("bla",3);  
  tm.setParam(param);
  
  TestMapping tm_copy;
  tm_copy = tm;
  
  TEST_EQUAL(tm_copy.getParam() == param,true)
RESULT

CHECK((BaseMapping(const BaseMapping& source)))
  TestMapping tm;
  Param param;
  param.setValue("bla",3);  
  tm.setParam(param);
  
  TestMapping tm_copy(tm);
  
  TEST_EQUAL(tm_copy.getParam() == param,true)
RESULT

CHECK((virtual const Param& getParam() const))
  TestMapping tm;
  Param param;
  
  TEST_EQUAL(tm.getParam() == param,true)
RESULT

CHECK((virtual const String getName()=0))
  
RESULT

CHECK((virtual void apply(DoubleReal &pos) const =0))
  
RESULT

CHECK((virtual void apply(DPosition< 1 > &) const =0))
  
RESULT

CHECK((virtual void setParam(const Param &p)))
  TestMapping tm;
  Param param;
  param.setValue("bla",3);  
  tm.setParam(param);
  
  TEST_EQUAL(tm.getParam() == param,true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



