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

#include <OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>

#include <vector>

///////////////////////////

START_TEST(LinearMapping, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


LinearMapping* d10_ptr = 0;
CHECK((LinearMapping()))
	d10_ptr = new LinearMapping;
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK((~LinearMapping()))
	delete d10_ptr;
RESULT

CHECK((LinearMapping(const LinearMapping& source)))

	LinearMapping c1(0.0, 10.0);
	LinearMapping c2(c1);
	
	TEST_EQUAL(c1.getSlope()==c2.getSlope(),true)	
	TEST_EQUAL(c1.getIntercept()==c2.getIntercept(),true)	
RESULT

CHECK((LinearMapping& operator = (const LinearMapping& source)))
	LinearMapping c1(0.0, 10.0);
	LinearMapping c2 = c1;
	
	TEST_EQUAL(c1.getSlope()==c2.getSlope(),true)	
	TEST_EQUAL(c1.getIntercept()==c2.getIntercept(),true)	
RESULT


CHECK((LinearMapping(DoubleReal slope, DoubleReal intercept)))
  LinearMapping c1(0.0, 10.0);
	LinearMapping c2 = c1;
	
	TEST_EQUAL(c1.getSlope()==c2.getSlope(),true)	
	TEST_EQUAL(c1.getIntercept()==c2.getIntercept(),true)	
RESULT

CHECK((void setParam(DoubleReal sl, DoubleReal in)))
  LinearMapping c1;
  c1.setParam(0.0, 10.0);
	
	TEST_REAL_EQUAL(c1.getSlope(),0)	
	TEST_REAL_EQUAL(c1.getIntercept(),10)	
RESULT

CHECK((void apply(DPosition< 1 > &pos) const))
  LinearMapping c1;
  c1.setParam(5, 2);
	DPosition<1> pos(1);
	c1.apply(pos);
	
	TEST_EQUAL(pos == 7,true)	
RESULT

CHECK((void apply(DoubleReal &pos) const))
  LinearMapping c1;
  c1.setParam(5, 2);
	DoubleReal pos(1);
	c1.apply(pos);
	
	TEST_REAL_EQUAL(pos,7)	
RESULT

CHECK((DoubleReal getSlope() const))
  LinearMapping c1;
  
	TEST_REAL_EQUAL(c1.getSlope(),1)	
RESULT

CHECK((void setSlope(DoubleReal sl)))
  LinearMapping c1;
  c1.setSlope(5);
	
	TEST_REAL_EQUAL(c1.getSlope(),5)	
RESULT

CHECK((DoubleReal getIntercept() const))
  LinearMapping c1;
  
	TEST_REAL_EQUAL(c1.getIntercept(),0)	
RESULT

CHECK((void setIntercept(DoubleReal in)))
  LinearMapping c1;
  c1.setIntercept(5);
	
	TEST_REAL_EQUAL(c1.getIntercept(),5)	
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
