// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

/////////////////////////////////////////////////////////////


#include <sstream>

using namespace OpenMS;
using namespace std;

START_TEST(DBoundingBox, "$Id$")

/////////////////////////////////////////////////////////////

typedef DBoundingBox<1> BB1;
typedef DBoundingBox<2> BB2;

BB1* ptr1 = nullptr;
BB1* nullPointer1 = nullptr;
START_SECTION(DBoundingBox())
	ptr1 = new BB1;
  TEST_NOT_EQUAL(ptr1, nullPointer1)
END_SECTION

START_SECTION(~DBoundingBox())
	delete ptr1;
END_SECTION

BB2* ptr2 = nullptr;
BB2* nullPointer2 = nullptr;
START_SECTION([EXTRA] DBoundingBox())
	ptr2 = new BB2;
  TEST_NOT_EQUAL(ptr2, nullPointer2)
END_SECTION

START_SECTION([EXTRA] ~DBoundingBox())
	delete ptr2;
END_SECTION

START_SECTION(DBoundingBox(const PositionType& minimum, const PositionType& maximum))
{
	DPosition<1> min(2), max(5);
	BB1 bb(min,max);
	TEST_REAL_SIMILAR(bb.minPosition()[0], 2);
	TEST_REAL_SIMILAR(bb.maxPosition()[0], 5);
}
END_SECTION

START_SECTION(DBoundingBox(const DBoundingBox &rhs))
	BB2 bb(DPosition<2>(1,2),DPosition<2>(3,4));
	BB2 bb_copy(bb);
	TEST_REAL_SIMILAR(bb.minPosition()[0],bb_copy.minPosition()[0]);
	TEST_REAL_SIMILAR(bb.minPosition()[1],bb_copy.minPosition()[1]);
	TEST_REAL_SIMILAR(bb.maxPosition()[0],bb_copy.maxPosition()[0]);
	TEST_REAL_SIMILAR(bb.maxPosition()[1],bb_copy.maxPosition()[1]);
END_SECTION

START_SECTION(DBoundingBox& operator=(const DBoundingBox &rhs))
	BB2 bb(DPosition<2>(1,2),DPosition<2>(3,4));
	BB2 bb_copy;
	bb_copy = bb;
	TEST_REAL_SIMILAR(bb.minPosition()[0],bb_copy.minPosition()[0]);
	TEST_REAL_SIMILAR(bb.minPosition()[1],bb_copy.minPosition()[1]);
	TEST_REAL_SIMILAR(bb.maxPosition()[0],bb_copy.maxPosition()[0]);
	TEST_REAL_SIMILAR(bb.maxPosition()[1],bb_copy.maxPosition()[1]);
END_SECTION

START_SECTION(DBoundingBox& operator=(const Base &rhs))
	BB2::Base bb(DPosition<2>(1,2),DPosition<2>(3,4));
	BB2 bb_copy;
	bb_copy = bb;
	TEST_REAL_SIMILAR(bb.minPosition()[0],bb_copy.minPosition()[0]);
	TEST_REAL_SIMILAR(bb.minPosition()[1],bb_copy.minPosition()[1]);
	TEST_REAL_SIMILAR(bb.maxPosition()[0],bb_copy.maxPosition()[0]);
	TEST_REAL_SIMILAR(bb.maxPosition()[1],bb_copy.maxPosition()[1]);
END_SECTION

START_SECTION(bool isEmpty() const)
	BB2 bb;
	bb = BB2::empty;
	TEST_EQUAL(bb.isEmpty(),true);
	bb = BB2::zero;
	TEST_EQUAL(bb.isEmpty(),true);
	bb = BB2(DPosition<2>(1,2),DPosition<2>(3,4));
	TEST_EQUAL(bb.isEmpty(),false);
END_SECTION

START_SECTION(void enlarge(const PositionType& p))
  BB2 bb2h;
  TEST_EQUAL(bb2h.encloses(11,13),false);
  TEST_EQUAL(bb2h.encloses(10,1),false);
  bb2h.enlarge(BB2::PositionType(11,13));
  TEST_EQUAL(bb2h.encloses(11,13),true);
  TEST_EQUAL(bb2h.encloses(10,1),false);
  bb2h.enlarge(BB2::PositionType(9,0));
  TEST_EQUAL(bb2h.encloses(11,13),true);
  TEST_EQUAL(bb2h.encloses(10,1),true);
END_SECTION

START_SECTION((void enlarge(CoordinateType x, CoordinateType y)))
  BB2 bb2h;
  TEST_EQUAL(bb2h.encloses(11,13),false);
  TEST_EQUAL(bb2h.encloses(10,1),false);
  bb2h.enlarge(11,13);
  TEST_EQUAL(bb2h.encloses(11,13),true);
  TEST_EQUAL(bb2h.encloses(10,1),false);
  bb2h.enlarge(9,0);
  TEST_EQUAL(bb2h.encloses(11,13),true);
  TEST_EQUAL(bb2h.encloses(10,1),true);
END_SECTION

START_SECTION(bool operator == (const DBoundingBox& rhs) const)
	BB2 bb2;
  bb2.enlarge(9,0);
	BB2 bb2_copy(bb2);
	TEST_EQUAL(bb2==bb2_copy,true);
END_SECTION

START_SECTION(bool operator == (const Base& rhs) const)
	BB2 bb2;
  bb2.enlarge(9,0);
	BB2::Base bb2_copy_base(bb2);
	TEST_EQUAL(bb2==bb2_copy_base,true);
END_SECTION

START_SECTION((bool encloses(CoordinateType x, CoordinateType y) const))
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

END_SECTION

START_SECTION(bool encloses(const PositionType& position) const)
BB2 tmp;
tmp.setMinX( 100 );
tmp.setMinY( 200 );
tmp.setMaxX( 300 );
tmp.setMaxY( 400 );
TEST_EQUAL(tmp.encloses(BB2::PositionType(10,200)),false);
TEST_EQUAL(tmp.encloses(BB2::PositionType(100,200)),true);
TEST_EQUAL(tmp.encloses(BB2::PositionType(200,200)),true);
TEST_EQUAL(tmp.encloses(BB2::PositionType(300,200)),true);
TEST_EQUAL(tmp.encloses(BB2::PositionType(310,200)),false);

TEST_EQUAL(tmp.encloses(BB2::PositionType(10,400)),false);
TEST_EQUAL(tmp.encloses(BB2::PositionType(100,400)),true);
TEST_EQUAL(tmp.encloses(BB2::PositionType(200,400)),true);
TEST_EQUAL(tmp.encloses(BB2::PositionType(300,400)),true);
TEST_EQUAL(tmp.encloses(BB2::PositionType(310,400)),false);

TEST_EQUAL(tmp.encloses(BB2::PositionType(200,190)),false);
TEST_EQUAL(tmp.encloses(BB2::PositionType(200,200)),true);
TEST_EQUAL(tmp.encloses(BB2::PositionType(200,300)),true);
TEST_EQUAL(tmp.encloses(BB2::PositionType(200,400)),true);
TEST_EQUAL(tmp.encloses(BB2::PositionType(200,410)),false);

TEST_EQUAL(tmp.encloses(BB2::PositionType(0,0)),false);

TEST_EQUAL(tmp.encloses(BB2::PositionType(100,200)),true);
TEST_EQUAL(tmp.encloses(BB2::PositionType(300,200)),true);
TEST_EQUAL(tmp.encloses(BB2::PositionType(100,400)),true);
TEST_EQUAL(tmp.encloses(BB2::PositionType(300,400)),true);
END_SECTION

START_SECTION(bool intersects(const DBoundingBox& bounding_box) const)
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
	r3.setMax(r2.maxPosition()+one);
	TEST_EQUAL(r2.intersects(r3),true)
	r3.setMin(r2.maxPosition()+one);
	r3.setMax(r2.maxPosition()+two);
	TEST_EQUAL(r2.intersects(r3),false)
	r3.setMin(r2.minPosition());
	r3.setMinX(10.0f);
	r3.setMax(r3.minPosition()+one);
	TEST_EQUAL(r2.intersects(r3),false)
	r3.setMinX(-10.0f);
	r3.setMinY(-10.0f);
	r3.setMax(r3.minPosition()+one);
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
	r3.setMax(r3.minPosition()+one);
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
END_SECTION

START_SECTION((template <UInt D> std::ostream & operator<<(std::ostream &os, const DBoundingBox< D > &bounding_box)))
{
	std::ostringstream os;
	DPosition<1> min(2), max(5);
	BB1 bb(min,max);
	os << bb;
  TEST_STRING_EQUAL( os.str(),
		"--DBOUNDINGBOX BEGIN--\n"
		"MIN --> 2.0\n"
		"MAX --> 5.0\n"
		"--DBOUNDINGBOX END--\n");
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


