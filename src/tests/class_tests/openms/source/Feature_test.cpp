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

#include <OpenMS/KERNEL/Feature.h>

///////////////////////////

START_TEST(Feature, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

typedef OpenMS::BaseFeature::QualityType QualityType;

Feature* d_ptr = nullptr;
Feature* nullPointer = nullptr;
START_SECTION((Feature()))
{
  d_ptr = new Feature;
  TEST_NOT_EQUAL(d_ptr, nullPointer);
}
END_SECTION

START_SECTION((~Feature()))
{
  delete d_ptr;
}
END_SECTION

START_SECTION((QualityType getOverallQuality() const))
  Feature p;
  TEST_REAL_SIMILAR(p.getOverallQuality(), 0.0)
  p.setOverallQuality((QualityType)123.456);
  TEST_REAL_SIMILAR(p.getOverallQuality(), 123.456)
  p.setOverallQuality((QualityType)-0.12345);
  TEST_REAL_SIMILAR(p.getOverallQuality(), -0.12345)
  p.setOverallQuality( (QualityType)0.0);
  TEST_REAL_SIMILAR(p.getOverallQuality(), 0.0)
END_SECTION

START_SECTION((void setOverallQuality(QualityType q)))
  Feature p;
  p.setOverallQuality((QualityType)123.456);
  TEST_REAL_SIMILAR(p.getOverallQuality(), 123.456)
  p.setOverallQuality( (QualityType)-0.12345);
  TEST_REAL_SIMILAR(p.getOverallQuality(), -0.12345)
  p.setOverallQuality( (QualityType)0.0);
  TEST_REAL_SIMILAR(p.getOverallQuality(), 0.0)
END_SECTION

START_SECTION((QualityType getQuality(Size index) const ))
  Feature p;
  TEST_REAL_SIMILAR(p.getQuality(0), 0.0)
  p.setQuality( 0, (QualityType)123.456);
  TEST_REAL_SIMILAR(p.getQuality(0), 123.456)
  p.setQuality( 0, (QualityType)-0.12345);
  TEST_REAL_SIMILAR(p.getQuality(0), -0.12345)
  p.setQuality( 0, (QualityType)0.0);
  TEST_REAL_SIMILAR(p.getQuality(0), 0.0)
  TEST_REAL_SIMILAR(p.getQuality(1), 0.0)
  TEST_PRECONDITION_VIOLATED(p.getQuality(10))
END_SECTION

START_SECTION((void setQuality(Size index, QualityType q)))
  Feature p;
  p.setQuality( 1, (QualityType)123.456);
  TEST_REAL_SIMILAR(p.getQuality(1), 123.456)
  p.setQuality( 1, (QualityType)-0.12345);
  TEST_REAL_SIMILAR(p.getQuality(1), -0.12345)
  p.setQuality( 1, (QualityType)0.0);
  TEST_REAL_SIMILAR(p.getQuality(0), 0.0)
  TEST_REAL_SIMILAR(p.getQuality(1), 0.0)
  TEST_PRECONDITION_VIOLATED(p.setQuality( 10, (QualityType)1.0))
END_SECTION

//do not change these datastructures, they are used in the following tests...
std::vector< ConvexHull2D > hulls(2);
hulls[0].addPoint(DPosition<2>(1.0,2.0));
hulls[0].addPoint(DPosition<2>(3.0,4.0));
hulls[1].addPoint(DPosition<2>(0.5,0.0));
hulls[1].addPoint(DPosition<2>(1.0,1.0));

START_SECTION((const vector<ConvexHull2D>& getConvexHulls() const))
  Feature tmp;
  TEST_EQUAL(tmp.getConvexHulls().size(),0)
END_SECTION

START_SECTION((vector<ConvexHull2D>& getConvexHulls()))
  Feature tmp;
  tmp.setConvexHulls(hulls);
  TEST_EQUAL(tmp.getConvexHulls().size(),2)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[0].getHullPoints()[0][0],1.0)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[0].getHullPoints()[0][1],2.0)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[0].getHullPoints()[1][0],3.0)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[0].getHullPoints()[1][1],4.0)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[1].getHullPoints()[0][0],0.5)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[1].getHullPoints()[0][1],0.0)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[1].getHullPoints()[1][0],1.0)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[1].getHullPoints()[1][1],1.0)
END_SECTION

START_SECTION((void setConvexHulls(const vector<ConvexHull2D>& hulls)))
  Feature tmp;
  tmp.setConvexHulls(hulls);
  TEST_EQUAL(tmp.getConvexHulls().size(),2)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[0].getHullPoints()[0][0],1.0)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[0].getHullPoints()[0][1],2.0)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[0].getHullPoints()[1][0],3.0)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[0].getHullPoints()[1][1],4.0)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[1].getHullPoints()[0][0],0.5)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[1].getHullPoints()[0][1],0.0)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[1].getHullPoints()[1][0],1.0)
  TEST_REAL_SIMILAR(tmp.getConvexHulls()[1].getHullPoints()[1][1],1.0)
END_SECTION

START_SECTION((ConvexHull2D& getConvexHull() const))
  Feature tmp;
  tmp.setConvexHulls(hulls);

  //check if the bounding box is ok
  DBoundingBox<2> bb = tmp.getConvexHull().getBoundingBox();
  TEST_REAL_SIMILAR(bb.minPosition()[0],0.5)
  TEST_REAL_SIMILAR(bb.minPosition()[1],0.0)
  TEST_REAL_SIMILAR(bb.maxPosition()[0],3.0)
  TEST_REAL_SIMILAR(bb.maxPosition()[1],4.0)

  //check the convex hull points
  TEST_EQUAL(tmp.getConvexHull().getHullPoints().size(),4)
  TEST_REAL_SIMILAR(tmp.getConvexHull().getHullPoints()[0][0],0.5)
  TEST_REAL_SIMILAR(tmp.getConvexHull().getHullPoints()[0][1],0.0)
  TEST_REAL_SIMILAR(tmp.getConvexHull().getHullPoints()[1][0],3.0)
  TEST_REAL_SIMILAR(tmp.getConvexHull().getHullPoints()[1][1],0.0)
  TEST_REAL_SIMILAR(tmp.getConvexHull().getHullPoints()[2][0],3.0)
  TEST_REAL_SIMILAR(tmp.getConvexHull().getHullPoints()[2][1],4.0)
  TEST_REAL_SIMILAR(tmp.getConvexHull().getHullPoints()[3][0],0.5)
  TEST_REAL_SIMILAR(tmp.getConvexHull().getHullPoints()[3][1],4.0)
END_SECTION

hulls[0].addPoint(DPosition<2>(3.0,2.0));
hulls[1].addPoint(DPosition<2>(2.0,1.0));

START_SECTION((bool encloses(double rt, double mz) const))
  Feature tmp;
  TEST_EQUAL(tmp.getConvexHull().getBoundingBox().isEmpty(), true)
  tmp.setConvexHulls(hulls);

  TEST_EQUAL(tmp.encloses(0.0,0.0), false);
  TEST_EQUAL(tmp.encloses(1.0,1.0), true);
  TEST_EQUAL(tmp.encloses(2.0,0.5), false);
  TEST_EQUAL(tmp.encloses(2.0,3.001), false);
  TEST_EQUAL(tmp.encloses(2.0,2.999), true);
  TEST_EQUAL(tmp.encloses(2.0,3.5), false);
  TEST_EQUAL(tmp.encloses(4.0,3.0), false);
  TEST_EQUAL(tmp.encloses(1.5,1.5), false);
  TEST_EQUAL(tmp.encloses(2.0,1.0), true);
  TEST_EQUAL(tmp.encloses(0.5,0.0), true);
  TEST_EQUAL(tmp.encloses(3.0,3.2), true);
END_SECTION

START_SECTION((Feature(const Feature &feature)))
{
  Feature::PositionType pos;
  pos[0] = 21.21;
  pos[1] = 22.22;
  Feature p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  p.setMetaValue("cluster_id",4711);
  p.setOverallQuality( (QualityType)0.9);
  p.setQuality( 0, (QualityType)0.1);
  p.setQuality( 1, (QualityType)0.2);
  p.setConvexHulls(hulls);
  p.getConvexHull(); //this pre-calculates the overall convex hull

  Feature::PositionType pos2;
  Feature::IntensityType i2;

  Feature copy_of_p(p);
  i2 = copy_of_p.getIntensity();
  pos2 = copy_of_p.getPosition();

  TEST_REAL_SIMILAR(i2, 123.456)

  TEST_REAL_SIMILAR(pos2[0], 21.21)
  TEST_REAL_SIMILAR(pos2[1], 22.22)

  TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"), DataValue(4711));

  Feature::QualityType q2;
  q2 = copy_of_p.getOverallQuality();
  TEST_REAL_SIMILAR(q2, 0.9)
  q2 = copy_of_p.getQuality(0);
  TEST_REAL_SIMILAR(q2, 0.1)
  q2 = copy_of_p.getQuality(1);
  TEST_REAL_SIMILAR(q2, 0.2)
  TEST_EQUAL(copy_of_p.getConvexHull().getHullPoints().size(),p.getConvexHull().getHullPoints().size())
  TEST_EQUAL(copy_of_p.getConvexHulls().size(),p.getConvexHulls().size())
}
END_SECTION

START_SECTION((Feature(const Feature&& source)))
{
  // Ensure that Feature has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(Feature(std::declval<Feature&&>())), true)

  Feature::PositionType pos;
  pos[0] = 21.21;
  pos[1] = 22.22;
  Feature p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  p.setMetaValue("cluster_id",4711);
  p.setOverallQuality( (QualityType)0.9);
  p.setQuality( 0, (QualityType)0.1);
  p.setQuality( 1, (QualityType)0.2);
  p.setConvexHulls(hulls);
  p.getConvexHull(); //this pre-calculates the overall convex hull

  Feature::PositionType pos2;
  Feature::IntensityType i2;

  //copy p so we can move one of them
  Feature orig = p;
  Feature copy_of_p(std::move(p));

  i2 = copy_of_p.getIntensity();
  pos2 = copy_of_p.getPosition();

  TEST_REAL_SIMILAR(i2, 123.456)

  TEST_REAL_SIMILAR(pos2[0], 21.21)
  TEST_REAL_SIMILAR(pos2[1], 22.22)

  TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"), DataValue(4711));

  Feature::QualityType q2;
  q2 = copy_of_p.getOverallQuality();
  TEST_REAL_SIMILAR(q2, 0.9)
  q2 = copy_of_p.getQuality(0);
  TEST_REAL_SIMILAR(q2, 0.1)
  q2 = copy_of_p.getQuality(1);
  TEST_REAL_SIMILAR(q2, 0.2)
  TEST_EQUAL(copy_of_p.getConvexHull().getHullPoints().size(), orig.getConvexHull().getHullPoints().size())
  TEST_EQUAL(copy_of_p.getConvexHulls().size(), orig.getConvexHulls().size())

  TEST_EQUAL(p.getConvexHull().getHullPoints().size(), 0)
  TEST_EQUAL(p.getConvexHulls().size(), 0)
}
END_SECTION

START_SECTION((Feature& operator = (const Feature& rhs)))
  Feature::PositionType pos;
  pos[0] = 21.21;
  pos[1] = 22.22;
  Feature p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  p.setOverallQuality( (QualityType)0.9);
  p.setQuality( 0, (QualityType)0.1);
  p.setQuality( 1, (QualityType)0.2);
  p.setMetaValue("cluster_id",4712);
  p.setConvexHulls(hulls);

  Feature::PositionType pos2;
  Feature::IntensityType i2;

  Feature copy_of_p;
  copy_of_p.getConvexHull(); //this pre-calculates the overall convex hull in order to check that the recalculation flag is copied correctly
  copy_of_p = p;

  i2 = copy_of_p.getIntensity();
  pos2 = copy_of_p.getPosition();

  Feature::QualityType q2;

  TEST_REAL_SIMILAR(i2, 123.456)
  TEST_REAL_SIMILAR(pos2[0], 21.21)
  TEST_REAL_SIMILAR(pos2[1], 22.22)
  q2 = copy_of_p.getOverallQuality();
  TEST_REAL_SIMILAR(q2, 0.9)
  q2 = copy_of_p.getQuality(0);
  TEST_REAL_SIMILAR(q2, 0.1)
  q2 = copy_of_p.getQuality(1);
  TEST_REAL_SIMILAR(q2, 0.2)
  TEST_EQUAL(copy_of_p.getConvexHull().getHullPoints().size(),p.getConvexHull().getHullPoints().size())
  TEST_EQUAL(copy_of_p.getConvexHulls().size(),p.getConvexHulls().size())
END_SECTION

START_SECTION((bool operator==(const Feature &rhs) const))
  Feature p1;
  Feature p2(p1);
  TEST_TRUE(p1 == p2)

  p1.setIntensity(5.0f);
  p1.setOverallQuality( (QualityType)0.9);
  p1.setQuality(0, (QualityType)0.1);
  TEST_EQUAL(p1==p2, false)
  p2.setIntensity(5.0f);
  p2.setOverallQuality( (QualityType)0.9);
  p2.setQuality(0, (QualityType)0.1);
  TEST_TRUE(p1 == p2)

  p1.getPosition()[0]=5;
  TEST_EQUAL(p1==p2, false)
  p2.getPosition()[0]=5;
  TEST_TRUE(p1 == p2)
END_SECTION

START_SECTION([EXTRA](Feature& operator != (const Feature& rhs)))
  Feature p1;
  Feature p2(p1);
  TEST_EQUAL(p1!=p2, false)

  p1.setIntensity(5.0f);
  TEST_FALSE(p1 == p2)
  p2.setIntensity(5.0f);
  TEST_EQUAL(p1!=p2, false)

  p1.getPosition()[0]=5;
  TEST_FALSE(p1 == p2)
  p2.getPosition()[0]=5;
  TEST_EQUAL(p1!=p2, false)
END_SECTION

START_SECTION(([EXTRA]meta info with copy constructor))
{
  Feature p;
  p.setMetaValue(2,String("bla"));
  Feature p2(p);
  TEST_EQUAL(p.getMetaValue(2), "bla")
  TEST_EQUAL(p2.getMetaValue(2), "bla")
  p.setMetaValue(2,String("bluff"));
  TEST_EQUAL(p.getMetaValue(2), "bluff")
  TEST_EQUAL(p2.getMetaValue(2), "bla")
}
END_SECTION

START_SECTION(([EXTRA]meta info with assignment))
{
  Feature p;
  p.setMetaValue(2,String("bla"));
  Feature p2 = p;
  TEST_EQUAL(p.getMetaValue(2), "bla")
  TEST_EQUAL(p2.getMetaValue(2), "bla")
  p.setMetaValue(2,String("bluff"));
  TEST_EQUAL(p.getMetaValue(2), "bluff")
  TEST_EQUAL(p2.getMetaValue(2), "bla")
}
END_SECTION

  
START_SECTION((std::vector<Feature>& getSubordinates()))
{
  // see below
  NOT_TESTABLE;
}
END_SECTION

START_SECTION((void setSubordinates(const std::vector<Feature>& rhs)))
{
  // see below
  NOT_TESTABLE;
}
END_SECTION

START_SECTION((const std::vector<Feature>& getSubordinates() const))
{
  Feature f1;
  f1.setRT(1001);
  f1.setMZ(1002);
  f1.setCharge(1003);
  Feature f1_cpy(f1);
  Feature f11;
  f11.setRT(1101);
  f11.setMZ(1102);
  Feature f12;
  f12.setRT(1201);
  f12.setMZ(1202);
  Feature f13;
  f13.setRT(1301);
  f13.setMZ(1302);
  TEST_EQUAL(f1.getSubordinates().empty(),true);
  f1.getSubordinates().push_back(f11);
  TEST_EQUAL(f1.getSubordinates().size(),1);
  f1.getSubordinates().push_back(f12);
  TEST_EQUAL(f1.getSubordinates().size(),2);
  f1.getSubordinates().push_back(f13);
  TEST_EQUAL(f1.getSubordinates().size(),3);
  TEST_EQUAL(f1.getRT(),1001);
  TEST_EQUAL(f1.getSubordinates()[0].getRT(),1101);
  TEST_EQUAL(f1.getSubordinates()[1].getRT(),1201);
  TEST_EQUAL(f1.getSubordinates()[2].getRT(),1301);
  const Feature &f1_cref = f1;
  TEST_EQUAL(f1_cref.getMZ(),1002);
  TEST_EQUAL(f1_cref.getSubordinates()[0].getMZ(),1102);
  TEST_EQUAL(f1_cref.getSubordinates()[1].getMZ(),1202);
  TEST_EQUAL(f1_cref.getSubordinates()[2].getMZ(),1302);
  TEST_NOT_EQUAL(f1_cref,f1_cpy);
  Feature f1_cpy2(f1);
  TEST_EQUAL(f1_cpy2,f1);
  f1.getSubordinates().clear();
  TEST_EQUAL(f1_cref,f1_cpy);

  Feature f2;
  f2.setRT(1001);
  f2.setMZ(1002);
  f2.setCharge(1003);
  TEST_NOT_EQUAL(f1_cpy2.getSubordinates().empty(),true);
  f2.setSubordinates(f1_cpy2.getSubordinates());
  TEST_EQUAL(f2,f1_cpy2);
}
END_SECTION

START_SECTION((template < typename Type > Size applyMemberFunction( Size (Type::*member_function)() )))
{
  Feature f;
  TEST_EQUAL(f.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId), 1);
  
  f.setUniqueId();
  TEST_EQUAL(f.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId), 0);
}
END_SECTION

START_SECTION((template < typename Type > Size applyMemberFunction( Size (Type::*member_function)() const) const))
{
  Feature f;
  TEST_EQUAL(f.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId), 1);
  
  f.setUniqueId();
  TEST_EQUAL(f.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId), 0);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
