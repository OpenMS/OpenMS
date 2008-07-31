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

// Includes go here....

#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

///////////////////////////

START_TEST(DBoundingBox, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

DBoundingBox<1>* ptr1 = 0;
CHECK(DBoundingBox())
	ptr1 = new DBoundingBox<1>;
	TEST_NOT_EQUAL(ptr1, 0)
RESULT

CHECK(~DBoundingBox())
	delete ptr1;
RESULT

DBoundingBox<2>* ptr2 = 0;
CHECK([EXTRA] DBoundingBox())
	ptr2 = new DBoundingBox<2>;
	TEST_NOT_EQUAL(ptr2, 0)
RESULT

CHECK([EXTRA] ~DBoundingBox())
	delete ptr2;
RESULT

// Misc stuff for testing
typedef DBoundingBox<1> BB1;
typedef DBoundingBox<2> BB2;
BB1 bb1;
BB2 bb2;

CHECK(bool isEmpty() const)
  TEST_EQUAL(bb1,BB1::empty)
  TEST_EQUAL(BB1::empty,bb1)
  TEST_EQUAL(bb2,BB2::empty)
RESULT

CHECK(void enlarge(const PositionType& p))
	bb2.enlarge(BB2::zero.min());
  TEST_EQUAL(bb2,BB2::zero);
	bb2.enlarge(BB2::empty.min());
  //STATUS(bb2);
	bb2.enlarge(BB2::empty.max());
  //STATUS(bb2);
RESULT

BB2 bb2h;

CHECK((void enlarge(CoordinateType x, CoordinateType y)))
  TEST_EQUAL(bb2h.encloses(11,13),false);
  TEST_EQUAL(bb2h.encloses(10,1),false);
  bb2h.enlarge(11,13);
  TEST_EQUAL(bb2h.encloses(11,13),true);
  TEST_EQUAL(bb2h.encloses(10,1),false);
  bb2h.enlarge(9,0);
  TEST_EQUAL(bb2h.encloses(11,13),true);
  TEST_EQUAL(bb2h.encloses(10,1),true);
RESULT

BB2 bb2f(bb2);

CHECK(bool operator == (const DBoundingBox& rhs) const)
	TEST_EQUAL(bb2f==bb2,true);
RESULT

CHECK(bool operator == (const Base& rhs) const)
	BB2::Base base(bb2);
	TEST_EQUAL(bb2f==base,true);
RESULT

CHECK(DBoundingBox& operator=(const DBoundingBox &rhs))
	bb2 = BB2::zero;
  TEST_EQUAL(bb2.isEmpty(),false);
  bb2 = BB2::empty;
  TEST_EQUAL(bb2.isEmpty(),true);
RESULT

CHECK(DBoundingBox& operator=(const Base &rhs))
	bb2 = BB2::Base::zero;
  TEST_EQUAL(bb2.isEmpty(),false);
  bb2 = BB2::Base::empty;
  TEST_EQUAL(bb2.isEmpty(),true);
RESULT

CHECK(DBoundingBox(const DBoundingBox& bounding_box))
	BB2 bb2b(bb2);
  TEST_EQUAL(bb2b,bb2);
RESULT

CHECK((bool encloses(CoordinateType x, CoordinateType y) const))
BB2 tmp;
tmp.setMinX( 100 );
tmp.setMinY( 200 );
tmp.setMaxX( 300 );
tmp.setMaxY( 400 );
TEST_EQUAL(tmp.encloses(10,200),false);
TEST_EQUAL(tmp.encloses(100,200),true);
TEST_EQUAL(tmp.encloses(200,200),true);
TEST_EQUAL(tmp.encloses(300,200),true);
TEST_EQUAL(tmp.encloses(310,200),false);

TEST_EQUAL(tmp.encloses(10,400),false);
TEST_EQUAL(tmp.encloses(100,400),true);
TEST_EQUAL(tmp.encloses(200,400),true);
TEST_EQUAL(tmp.encloses(300,400),true);
TEST_EQUAL(tmp.encloses(310,400),false);

TEST_EQUAL(tmp.encloses(200,190),false);
TEST_EQUAL(tmp.encloses(200,200),true);
TEST_EQUAL(tmp.encloses(200,300),true);
TEST_EQUAL(tmp.encloses(200,400),true);
TEST_EQUAL(tmp.encloses(200,410),false);

TEST_EQUAL(tmp.encloses(0,0),false);

TEST_EQUAL(tmp.encloses(100,200),true);
TEST_EQUAL(tmp.encloses(300,200),true);
TEST_EQUAL(tmp.encloses(100,400),true);
TEST_EQUAL(tmp.encloses(300,400),true);

RESULT

CHECK(bool encloses(const PositionType& position) const)
  // see above :-P
RESULT

CHECK(bool intersects(const DBoundingBox& bounding_box) const)
	DPosition<2> p1,p2,p3,one,two;
	p1[0]=-1.0f;
	p1[1]=-2.0f;
	p2[0]=3.0f;
	p2[1]=4.0f;
	p3[0]=-10.0f;
	p3[1]=20.0f;
	one[0]=1;
	one[1]=1;
	two[0]=2;
	two[1]=2;

	BB2 r2;
	r2.setMin(p1);
	r2.setMax(p2);
	BB2 r3(r2);
	TEST_EQUAL(r2.intersects(r3),true)
	r3.setMaxX(10.0f);
	TEST_EQUAL(r2.intersects(r3),true)
	r3.setMax(r2.max()+one);
	TEST_EQUAL(r2.intersects(r3),true)
	r3.setMin(r2.max()+one);
	r3.setMax(r2.max()+two);
	TEST_EQUAL(r2.intersects(r3),false)
	r3.setMin(r2.min());
	r3.setMinX(10.0f);
	r3.setMax(r3.min()+one);
	TEST_EQUAL(r2.intersects(r3),false)
	r3.setMinX(-10.0f);
	r3.setMinY(-10.0f);
	r3.setMax(r3.min()+one);
	TEST_EQUAL(r2.intersects(r3),false)
	r3.setMinX(-10.0f);
	r3.setMinY(-10.0f);
	r3.setMaxX(0.0f);
	r3.setMaxY(-9.0f);
	TEST_EQUAL(r2.intersects(r3),false)
	r3.setMinX(-10.0f);
	r3.setMinY(-10.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(-9.0f);
	TEST_EQUAL(r2.intersects(r3),false)
	r3.setMinX(-10.0f);
	r3.setMinY(0.0f);
	r3.setMaxX(-9.0f);
	r3.setMaxY(1.0f);
	TEST_EQUAL(r2.intersects(r3),false)
	r3.setMinX(-10.0f);
	r3.setMinY(10.0f);
	r3.setMax(r3.min()+one);
	TEST_EQUAL(r2.intersects(r3),false)
	r3.setMinX(-10.0f);
	r3.setMinY(0.0f);
	r3.setMaxX(-9.0f);
	r3.setMaxY(10.0f);
	TEST_EQUAL(r2.intersects(r3),false)
	r3.setMinX(9.0f);
	r3.setMinY(0.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(10.0f);
	TEST_EQUAL(r2.intersects(r3),false)
	r3.setMinX(9.0f);
	r3.setMinY(0.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(10.0f);
	TEST_EQUAL(r2.intersects(r3),false)
	r3.setMinX(9.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.intersects(r3),false)
	r3.setMinX(9.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(5.0f);
	TEST_EQUAL(r2.intersects(r3),false)
	r3.setMinX(-5.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(0.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.intersects(r3),true)
	r3.setMinX(-5.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(5.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.intersects(r3),true)
	r3.setMinX(-5.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(5.0f);
	r3.setMaxY(5.0f);
	TEST_EQUAL(r2.intersects(r3),true)
	r3.setMinX(0.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(0.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.intersects(r3),true)
	r3.setMinX(0.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(5.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.intersects(r3),true)
	r3.setMinX(0.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(5.0f);
	r3.setMaxY(5.0f);
	TEST_EQUAL(r2.intersects(r3),true)
RESULT

CHECK(DBoundingBox(const PositionType& minimum, const PositionType& maximum))
{
	DPosition<1> min(2), max(5);
	BB1 bb(min,max);
	TEST_REAL_EQUAL(bb.min()[0], 2);
	TEST_REAL_EQUAL(bb.max()[0], 5);
}
RESULT

CHECK(DBoundingBox(const DBoundingBox &rhs))
{
	DPosition<1> min(2), max(5);
	BB1 bb(min,max);
	BB1 bb_copy(bb);
	TEST_REAL_EQUAL(bb_copy.min()[0], 2);
	TEST_REAL_EQUAL(bb_copy.max()[0], 5);
}
RESULT

CHECK((template <UInt D> std::ostream & operator<<(std::ostream &os, const DBoundingBox< D > &bounding_box)))
{
	std::ostringstream os;
	DPosition<1> min(2), max(5);
	BB1 bb(min,max);
	os << bb;
  TEST_STRING_EQUAL( os.str(),
		"--DBOUNDINGBOX BEGIN--\n"
		"MIN --> 2\n"
		"MAX --> 5\n"
		"--DBOUNDINGBOX END--\n");
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


