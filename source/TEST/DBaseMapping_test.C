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
#include <OpenMS/ANALYSIS/MAPMATCHING/DBaseMapping.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class TestMapping : public DBaseMapping<2>
{
  public:
  TestMapping() : DBaseMapping<2>(){}
  TestMapping(const TestMapping& tm) : DBaseMapping<2>(tm){}
  TestMapping& operator=(const TestMapping& tm)
  {
     DBaseMapping<2>::operator=(tm);
     return *this;
  }
    
   virtual void apply(DPosition<2>& ) const {}

   virtual void apply(KernelTraits::RealType& ) const {};
  
   virtual const String getName() { return "";}
};

START_TEST(DBaseMapping, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestMapping* ptr = 0;
CHECK(DBaseMapping())
	ptr = new TestMapping();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~DBaseMapping())
	delete ptr;
RESULT

CHECK(DBaseMapping& operator = (const DBaseMapping& rhs))
  TestMapping tm;
  Param param;
  param.setValue("bla",3);  
  tm.setParam(param);
  
  TestMapping tm_copy;
  tm_copy = tm;
  
  TEST_EQUAL(tm_copy.getParam() == param,true)
RESULT

CHECK(DBaseMapping(const DBaseMapping& source))
  TestMapping tm;
  Param param;
  param.setValue("bla",3);  
  tm.setParam(param);
  
  TestMapping tm_copy(tm);
  
  TEST_EQUAL(tm_copy.getParam() == param,true)
RESULT

CHECK(bool operator != (const DBaseMapping& rhs))
  TestMapping tm;
  Param param;
  param.setValue("bla",3);  
  
  TestMapping tm2;
  tm2.setParam(param);
  
  TEST_EQUAL(tm != tm2, true)
RESULT

CHECK(bool operator == (const DBaseMapping& rhs))
  TestMapping tm;
  Param param;
  param.setValue("bla",3);  
  tm.setParam(param);
  
  TestMapping tm2;
  tm2.setParam(param);
  
  TEST_EQUAL(tm == tm2, true)
RESULT

CHECK(const Param& getParam() const)
  TestMapping tm;
  Param param;
  
  TEST_EQUAL(tm.getParam() == param,true)
RESULT

CHECK(const String getName())
  
RESULT

CHECK(void apply( typename Traits::RealType & pos) const)
  
RESULT

CHECK(void apply(DPosition<D>& ) const)
  
RESULT

CHECK(void setParam(const Param& p))
  TestMapping tm;
  Param param;
  param.setValue("bla",3);  
  tm.setParam(param);
  
  TEST_EQUAL(tm.getParam() == param,true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



