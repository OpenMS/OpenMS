// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>

///////////////////////////

START_TEST(ConvexHull2D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ConvexHull2D* ptr = nullptr;
ConvexHull2D* nullPointer = nullptr;
START_SECTION((ConvexHull2D()))
	ptr = new ConvexHull2D;
	TEST_NOT_EQUAL(ptr, nullPointer)
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

START_SECTION((ConvexHull2D(const ConvexHull2D&& source)))
{
#ifndef OPENMS_COMPILER_MSVC
  // Ensure that ConvexHull2D has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  // Note that MSVS does not support noexcept move constructors for STL
  // constructs such as std::map.
  TEST_EQUAL(noexcept(ConvexHull2D(std::declval<ConvexHull2D&&>())), true)
#endif
}
END_SECTION

START_SECTION((void addPoints(const PointArrayType &points)))
	ConvexHull2D tmp;
	TEST_EQUAL(tmp.getHullPoints().size(),0)
	tmp.addPoints(vec);
	TEST_EQUAL(!tmp.getHullPoints().empty(),true)
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

START_SECTION((Size compress()))
{
  ConvexHull2D tmp;

  tmp.addPoint(DPosition<2>(1.,1.));
  tmp.addPoint(DPosition<2>(1.,10.));

  tmp.addPoint(DPosition<2>(2.,1.));
  tmp.addPoint(DPosition<2>(2.,10.));

  tmp.addPoint(DPosition<2>(3.,1.));
  tmp.addPoint(DPosition<2>(3.,10.));

  DBoundingBox<2> beforeCompress = tmp.getBoundingBox();

  TEST_EQUAL(tmp.compress() , 1)

  // second call should remove no points
  TEST_EQUAL(tmp.compress() , 0)


  TEST_EQUAL(tmp.getBoundingBox(), beforeCompress)

  tmp.addPoint(DPosition<2>(4.,1.));
  tmp.addPoint(DPosition<2>(4.,10.));

  tmp.addPoint(DPosition<2>(5.,2.));
  tmp.addPoint(DPosition<2>(5.,10.));

  tmp.addPoint(DPosition<2>(6.,1.));
  tmp.addPoint(DPosition<2>(6.,10.));

  beforeCompress = tmp.getBoundingBox();

  TEST_EQUAL(tmp.compress() , 1)

  // second call should remove no points
  TEST_EQUAL(tmp.compress() , 0)

  TEST_EQUAL(tmp.getBoundingBox(), beforeCompress)

  // check if encloses still works correct

  TEST_EQUAL(tmp.encloses(DPosition<2>(1.1, 5.)), true)
  TEST_EQUAL(tmp.encloses(DPosition<2>(2.1, 5.)), true)
  TEST_EQUAL(tmp.encloses(DPosition<2>(3.1, 5.)), true)
  TEST_EQUAL(tmp.encloses(DPosition<2>(4.1, 5.)), true)
  TEST_EQUAL(tmp.encloses(DPosition<2>(5.1, 5.)), true)
  TEST_EQUAL(tmp.encloses(DPosition<2>(5.1, 1.)), false)
  TEST_EQUAL(tmp.encloses(DPosition<2>(5.9, 5.)), true)
}
END_SECTION


START_SECTION((void expandToBoundingBox()))
{
  ConvexHull2D tmp;

  tmp.addPoint(DPosition<2>(1.,1.));
  tmp.addPoint(DPosition<2>(1.,10.));
  tmp.addPoint(DPosition<2>(2.,1.));
  tmp.addPoint(DPosition<2>(2.,10.));
  tmp.addPoint(DPosition<2>(3.,1.));
  tmp.addPoint(DPosition<2>(3.,10.));
  tmp.addPoint(DPosition<2>(4.,1.));
  tmp.addPoint(DPosition<2>(4.,10.));
  tmp.addPoint(DPosition<2>(5.,2.));
  tmp.addPoint(DPosition<2>(5.,10.));
  tmp.addPoint(DPosition<2>(6.,1.));
  tmp.addPoint(DPosition<2>(6.,10.));

  ConvexHull2D original(tmp);

	// Make sure we are left with only four points afterwards.
	tmp.expandToBoundingBox();
	TEST_EQUAL(tmp.getHullPoints().size(), 4)

  // second call should remove no points
	tmp.expandToBoundingBox();
	TEST_EQUAL(tmp.getHullPoints().size(), 4)

	// Check that values agree with min/max of the 
	// enclosed points.
	float min_x, min_y, max_x, max_y;
	min_x = tmp.getHullPoints()[0][0];
	min_y = tmp.getHullPoints()[0][1];
	max_x = min_x;
	max_y = min_y;
	for (Size i = 0; i < tmp.getHullPoints().size(); ++i)
	{
		float x = tmp.getHullPoints()[i][0];
		float y = tmp.getHullPoints()[i][1];
		min_x = std::min(min_x, x);
		max_x = std::max(max_x, x);
		min_y = std::min(min_y, y);
		max_y = std::max(max_y, y);
	}
	float o_min_x, o_min_y, o_max_x, o_max_y;
	o_min_x = original.getHullPoints()[0][0];
	o_min_y = original.getHullPoints()[0][1];
	o_max_x = o_min_x;
	o_max_y = o_min_y;
	for (Size i = 0; i < original.getHullPoints().size(); ++i)
	{
		float x = original.getHullPoints()[i][0];
		float y = original.getHullPoints()[i][1];
		o_min_x = std::min(o_min_x, x);
		o_max_x = std::max(o_max_x, x);
		o_min_y = std::min(o_min_y, y);
		o_max_y = std::max(o_max_y, y);
	}
	TEST_REAL_SIMILAR(min_x, o_min_x)
	TEST_REAL_SIMILAR(min_y, o_min_y)
	TEST_REAL_SIMILAR(max_x, o_max_x)
	TEST_REAL_SIMILAR(max_y, o_max_y)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
