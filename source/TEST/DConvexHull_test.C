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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DConvexHull.h>
///////////////////////////

START_TEST(DRange<D>, "$id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

DConvexHull<2>* ptr = 0;
CHECK((DConvexHull()))
	ptr = new DConvexHull<2>;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(([EXTRA] ~DConvexHull()))
	delete ptr;
RESULT

CHECK((const PointArrayType& getPoints() const))
	DConvexHull<2> tmp;
	TEST_EQUAL(tmp.getPoints().size(),0)
RESULT

//do not change these definitions, they are used in many tests
DPosition<2> p1(1.0,2.0);
DPosition<2> p2(3.0,4.0);
DPosition<2> p3(5.0,0.0);

DPosition<2> p4(1.0,1.0);
DPosition<2> p5(3.0,1.0);
DPosition<2> p6(1.0,3.0);

vector<DPosition<2> > vec;
vec.push_back(p1);
vec.push_back(p2);
vec.push_back(p3);

vector<DPosition<2> > vec2;
vec2.push_back(p4);
vec2.push_back(p5);
vec2.push_back(p6);

CHECK((DConvexHull& operator=(const PointArrayType& points)))
	DConvexHull<2> tmp;
	vector<DPosition<2> > vec3;
	vec3.push_back(p1);
	tmp = vec3;
	TEST_EQUAL(tmp.getPoints().size(),1)
	vec3.push_back(p2);
	tmp = vec3;
	TEST_EQUAL(tmp.getPoints().size(),2)
	vec3.push_back(p3);
	tmp = vec3;
	TEST_EQUAL(tmp.getPoints().size(),3)
RESULT

CHECK((DConvexHull& operator=(const DConvexHull& rhs)))
	DConvexHull<2> tmp,tmp2;
	tmp = vec;
	tmp2 = tmp;
	TEST_EQUAL(tmp2.getPoints().size(),3)
RESULT

CHECK((void clear()))
	DConvexHull<2> tmp;
	tmp = vec;
	tmp.clear();
	TEST_EQUAL(tmp.getPoints().size(),0)	
RESULT

CHECK((bool encloses(const PointType& point) const))
	DConvexHull<2> tmp;
	tmp=vec2;
	TEST_EQUAL(tmp.encloses(DPosition<2>(3.0,3.0)),false)
	TEST_EQUAL(tmp.encloses(DPosition<2>(0.0,0.0)),false)
	TEST_EQUAL(tmp.encloses(DPosition<2>(6.0,0.0)),false)
	TEST_EQUAL(tmp.encloses(DPosition<2>(0.0,6.0)),false)
	TEST_EQUAL(tmp.encloses(DPosition<2>(1.5,1.5)),true)
	TEST_EQUAL(tmp.encloses(DPosition<2>(1.0,1.0)),true)
	TEST_EQUAL(tmp.encloses(DPosition<2>(1.1,1.0)),true)
	TEST_EQUAL(tmp.encloses(DPosition<2>(1.2,2.5)),true)
	TEST_EQUAL(tmp.encloses(DPosition<2>(2.5,1.2)),true)
RESULT

CHECK((bool operator==(const DConvexHull& rhs) const))
	DConvexHull<2> tmp,tmp2;
	tmp=vec2;
	TEST_EQUAL(tmp==tmp2,false)
	tmp2=vec;
	TEST_EQUAL(tmp==tmp2,false)
	tmp2=vec2;
	TEST_EQUAL(tmp==tmp2,true)
RESULT

CHECK((DBoundingBox<D> getBoundingBox() const))
	//empty
	DConvexHull<2> tmp2;
	tmp2 = vec;
	DBoundingBox<2> bb2 = tmp2.getBoundingBox();
	TEST_REAL_EQUAL(bb2.min()[0],1.0)
	TEST_REAL_EQUAL(bb2.min()[1],0.0)
	TEST_REAL_EQUAL(bb2.max()[0],5.0)
	TEST_REAL_EQUAL(bb2.max()[1],4.0)
	
	//full
	DConvexHull<2> tmp;
	DBoundingBox<2> bb;

	bb = tmp.getBoundingBox();
	TEST_REAL_EQUAL(bb.isEmpty(),true)
	
	tmp = vec2;
	bb = tmp.getBoundingBox();
	TEST_REAL_EQUAL(bb.min()[0],1.0)
	TEST_REAL_EQUAL(bb.min()[1],1.0)
	TEST_REAL_EQUAL(bb.max()[0],3.0)
	TEST_REAL_EQUAL(bb.max()[1],3.0)

	tmp = vec;
	bb = tmp.getBoundingBox();
	TEST_REAL_EQUAL(bb.min()[0],1.0)
	TEST_REAL_EQUAL(bb.min()[1],0.0)
	TEST_REAL_EQUAL(bb.max()[0],5.0)
	TEST_REAL_EQUAL(bb.max()[1],4.0)

	vector<DPosition<2> > vec3;
	vec3.push_back(p1);
	tmp = vec3;
	bb = tmp.getBoundingBox();
	TEST_REAL_EQUAL(bb.min()[0],1.0)
	TEST_REAL_EQUAL(bb.min()[1],2.0)
	TEST_REAL_EQUAL(bb.max()[0],1.0)
	TEST_REAL_EQUAL(bb.max()[1],2.0)

	vec3.push_back(p2);
	tmp = vec3;
	bb = tmp.getBoundingBox();
	TEST_REAL_EQUAL(bb.min()[0],1.0)
	TEST_REAL_EQUAL(bb.min()[1],2.0)
	TEST_REAL_EQUAL(bb.max()[0],3.0)
	TEST_REAL_EQUAL(bb.max()[1],4.0)
RESULT

CHECK((bool addPoint(const PointType& point)))
	DConvexHull<2> tmp;
	tmp=vec2;
	TEST_EQUAL(tmp.addPoint(DPosition<2>(1.5,1.5)),false)
	TEST_EQUAL(tmp.addPoint(DPosition<2>(1.0,1.0)),false)
	TEST_EQUAL(tmp.addPoint(DPosition<2>(3.0,2.5)),true)
	TEST_EQUAL(tmp.addPoint(DPosition<2>(0.5,0.5)),true)	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
