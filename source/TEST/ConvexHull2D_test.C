// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>

///////////////////////////

START_TEST(DRange<D>, "$id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ConvexHull2D* ptr = 0;
START_SECTION((ConvexHull2D()))
	ptr = new ConvexHull2D;
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(([EXTRA] ~ConvexHull2D()))
	delete ptr;
END_SECTION

START_SECTION((const PointArrayType& getHullPoints() const))
	ConvexHull2D tmp;
	TEST_EQUAL(tmp.getHullPoints().size(),0)
END_SECTION

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

START_SECTION(void setHullPoints(const PointArrayType& points))
	ConvexHull2D tmp;
	vector<DPosition<2> > vec3;
	vec3.push_back(p1);
	tmp.setHullPoints(vec3);
	TEST_EQUAL(tmp.getHullPoints().size(),1)

	vec3.push_back(p2);
	tmp.setHullPoints(vec3);
	TEST_EQUAL(tmp.getHullPoints().size(),2)

	vec3.push_back(p3);
	tmp.setHullPoints(vec3);
	TEST_EQUAL(tmp.getHullPoints().size(),3)

	vec3.push_back(p5);
	tmp.setHullPoints(vec3);
	TEST_EQUAL(tmp.getHullPoints().size(),4)	
END_SECTION

START_SECTION((ConvexHull2D& operator=(const ConvexHull2D& rhs)))
	ConvexHull2D tmp,tmp2;
	tmp.setHullPoints(vec);
	tmp2 = tmp;
	TEST_EQUAL(tmp2.getHullPoints().size(),3)
END_SECTION


START_SECTION((void addPoints(const PointArrayType &points)))
	ConvexHull2D tmp;
	TEST_EQUAL(tmp.getHullPoints().size(),0)
	tmp.addPoints(vec);
	TEST_EQUAL(tmp.getHullPoints().size()!=0,true)
END_SECTION



START_SECTION((void clear()))
	vector<DPosition<2> > vec3;
	vec3.push_back(p1);
	vec3.push_back(p2);
	vec3.push_back(p3);
	vec3.push_back(p5);
	ConvexHull2D tmp;
	tmp.setHullPoints(vec3);
	TEST_EQUAL(tmp.getHullPoints().size(),4)
	tmp.clear();
	TEST_EQUAL(tmp.getHullPoints().size(),0)	
	
	tmp.addPoints(vec3);
	TEST_EQUAL(tmp.getHullPoints().size(),4)
	tmp.clear();
	TEST_EQUAL(tmp.getHullPoints().size(),0)	

END_SECTION

START_SECTION((bool encloses(const PointType& point) const))
	ConvexHull2D tmp;
	// setting hull points alone does not allow to query encloses()
	tmp.setHullPoints(vec2);
	TEST_EXCEPTION(Exception::NotImplemented, tmp.encloses(DPosition<2>(1.0,1.0)))

	tmp.addPoints(vec);
	tmp.addPoints(vec2);
	TEST_EQUAL(tmp.encloses(DPosition<2>(3.0,3.0)),true)
	TEST_EQUAL(tmp.encloses(DPosition<2>(0.0,0.0)),false)
	TEST_EQUAL(tmp.encloses(DPosition<2>(6.0,0.0)),false)
	TEST_EQUAL(tmp.encloses(DPosition<2>(0.0,6.0)),false)
	TEST_EQUAL(tmp.encloses(DPosition<2>(1.5,1.5)),true)
	TEST_EQUAL(tmp.encloses(DPosition<2>(1.0,1.0)),true)
	TEST_EQUAL(tmp.encloses(DPosition<2>(1.1,1.0)),true)
	TEST_EQUAL(tmp.encloses(DPosition<2>(1.2,2.5)),true)
	TEST_EQUAL(tmp.encloses(DPosition<2>(1.2,3.21)),false)
	TEST_EQUAL(tmp.encloses(DPosition<2>(1.4,0.99)),false)
	TEST_EQUAL(tmp.encloses(DPosition<2>(2.5,1.2)),true)
	TEST_EQUAL(tmp.encloses(DPosition<2>(1.0,1.1)),true)
	TEST_EQUAL(tmp.encloses(DPosition<2>(3.0,1.0)),true)
	TEST_EQUAL(tmp.encloses(DPosition<2>(5.0,0.0)),true)
END_SECTION

START_SECTION((bool operator==(const ConvexHull2D& rhs) const))
	ConvexHull2D tmp,tmp2;
	tmp.setHullPoints(vec2);
	TEST_EQUAL(tmp==tmp2,false)
	tmp2.setHullPoints(vec);
	TEST_EQUAL(tmp==tmp2,false)
	tmp2.setHullPoints(vec2);
	TEST_EQUAL(tmp==tmp2,true)
	tmp2.addPoints(vec);
	TEST_EQUAL(tmp==tmp2,false)
	tmp.addPoints(vec);
	TEST_EQUAL(tmp==tmp2,true)
END_SECTION

START_SECTION((DBoundingBox<2> getBoundingBox() const))
	//empty
	ConvexHull2D tmp2;
	TEST_EQUAL(tmp2.getBoundingBox().isEmpty(), true);
	tmp2.setHullPoints(vec);
	DBoundingBox<2> bb2 = tmp2.getBoundingBox();
	TEST_REAL_SIMILAR(bb2.minPosition()[0],1.0)
	TEST_REAL_SIMILAR(bb2.minPosition()[1],0.0)
	TEST_REAL_SIMILAR(bb2.maxPosition()[0],5.0)
	TEST_REAL_SIMILAR(bb2.maxPosition()[1],4.0)
	
	//full
	ConvexHull2D tmp;
	DBoundingBox<2> bb;

	bb = tmp.getBoundingBox();
	TEST_EQUAL(bb.isEmpty(),true)
	
	tmp.setHullPoints(vec2);
	bb = tmp.getBoundingBox();
	TEST_REAL_SIMILAR(bb.minPosition()[0],1.0)
	TEST_REAL_SIMILAR(bb.minPosition()[1],1.0)
	TEST_REAL_SIMILAR(bb.maxPosition()[0],3.0)
	TEST_REAL_SIMILAR(bb.maxPosition()[1],3.0)

	tmp.setHullPoints(vec);
	bb = tmp.getBoundingBox();
	TEST_REAL_SIMILAR(bb.minPosition()[0],1.0)
	TEST_REAL_SIMILAR(bb.minPosition()[1],0.0)
	TEST_REAL_SIMILAR(bb.maxPosition()[0],5.0)
	TEST_REAL_SIMILAR(bb.maxPosition()[1],4.0)

	vector<DPosition<2> > vec3;
	vec3.push_back(p1);
	tmp.setHullPoints(vec3);
	bb = tmp.getBoundingBox();
	TEST_REAL_SIMILAR(bb.minPosition()[0],1.0)
	TEST_REAL_SIMILAR(bb.minPosition()[1],2.0)
	TEST_REAL_SIMILAR(bb.maxPosition()[0],1.0)
	TEST_REAL_SIMILAR(bb.maxPosition()[1],2.0)

	vec3.push_back(p2);
	tmp.setHullPoints(vec3);
	bb = tmp.getBoundingBox();
	TEST_REAL_SIMILAR(bb.minPosition()[0],1.0)
	TEST_REAL_SIMILAR(bb.minPosition()[1],2.0)
	TEST_REAL_SIMILAR(bb.maxPosition()[0],3.0)
	TEST_REAL_SIMILAR(bb.maxPosition()[1],4.0)
END_SECTION

START_SECTION((bool addPoint(const PointType& point)))
	ConvexHull2D tmp;
	TEST_EQUAL(tmp.addPoint(DPosition<2>(1.5,1.5)),true)
	TEST_EQUAL(tmp.addPoint(DPosition<2>(1.0,1.0)),true)
	TEST_EQUAL(tmp.addPoint(DPosition<2>(1.0,1.5)),true)
	TEST_EQUAL(tmp.addPoint(DPosition<2>(1.0,1.2)),false)
	TEST_EQUAL(tmp.addPoint(DPosition<2>(3.0,2.5)),true)
	TEST_EQUAL(tmp.addPoint(DPosition<2>(3.0,1.5)),true)
	TEST_EQUAL(tmp.addPoint(DPosition<2>(3.0,2.5)),false)
	TEST_EQUAL(tmp.addPoint(DPosition<2>(3.0,2.0)),false)
	TEST_EQUAL(tmp.addPoint(DPosition<2>(0.5,0.5)),true)	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
