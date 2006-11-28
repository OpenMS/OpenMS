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
CHECK(DConvexHull())
	ptr = new DConvexHull<2>;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~DConvexHull())
	delete ptr;
RESULT

CHECK(const PointType& getPoints() const)
	DConvexHull<2> tmp;
	TEST_EQUAL(tmp.getPoints().size(),0)
RESULT

//do not change these points, they are used in many tests
DPosition<2> p1;
p1[0] = 1.0;
p1[1] = 2.0;
DPosition<2> p2;
p2[0] = 3.0;
p2[1] = 4.0;

CHECK(void addPoint(cost PointType& point))
	DConvexHull<2> tmp;
	tmp.addPoint(p1);
	TEST_EQUAL(tmp.getPoints().size(),1)
	TEST_REAL_EQUAL(tmp.getPoints()[0][0],1.0)
	TEST_REAL_EQUAL(tmp.getPoints()[0][1],2.0)
RESULT

CHECK(void clear())
	DConvexHull<2> tmp;
	tmp.addPoint(p1);
	tmp.clear();
	TEST_EQUAL(tmp.getPoints().size(),0)	
RESULT

CHECK(ConvexHull<D>& operator=(const ConvexHull& rhs))
	DConvexHull<2> tmp,tmp2;
	tmp.addPoint(p1);
	tmp2 = tmp;
	TEST_EQUAL(tmp2.getPoints().size(),1)
	TEST_REAL_EQUAL(tmp2.getPoints()[0][0],1.0)
	TEST_REAL_EQUAL(tmp2.getPoints()[0][1],2.0)
RESULT

CHECK(ConvexHull<D>& operator=(const PointArrayType& points))
	DConvexHull<2> tmp;
	vector<DPosition<2> > vec;
	vec.push_back(p1);
	vec.push_back(p2);
	tmp = vec;
	TEST_EQUAL(tmp.getPoints().size(),2)
	TEST_REAL_EQUAL(tmp.getPoints()[0][0],1.0)
	TEST_REAL_EQUAL(tmp.getPoints()[0][1],2.0)
	TEST_REAL_EQUAL(tmp.getPoints()[1][0],3.0)
	TEST_REAL_EQUAL(tmp.getPoints()[1][1],4.0)
RESULT

CHECK(DBoundingBox<D> getBoundingBox() const)
	DConvexHull<2> tmp;
	tmp.addPoint(p2);
	tmp.addPoint(p1);
	tmp.addPoint(p2);
	tmp.addPoint(p1);
	DBoundingBox<2> bb = tmp.getBoundingBox();
	TEST_REAL_EQUAL(bb.min()[0],1.0)
	TEST_REAL_EQUAL(bb.min()[1],2.0)
	TEST_REAL_EQUAL(bb.max()[0],3.0)
	TEST_REAL_EQUAL(bb.max()[1],4.0)
RESULT

CHECK(bool encloses(const PointType& point) const)
	//TODO
RESULT

CHECK(bool operator==(const DConvexHull& rhs) const)
	//TODO
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
