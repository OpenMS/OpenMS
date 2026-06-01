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

#include <OpenMS/KERNEL/RichPeak2D.h>

///////////////////////////

START_TEST(RichPeak2D<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

RichPeak2D* d10_ptr = nullptr;
RichPeak2D* d10_nullPointer = nullptr;
START_SECTION((RichPeak2D()))
  d10_ptr = new RichPeak2D;
  TEST_NOT_EQUAL(d10_ptr, d10_nullPointer)
END_SECTION

START_SECTION((~RichPeak2D()))
  delete d10_ptr;
END_SECTION

START_SECTION((RichPeak2D(const RichPeak2D &p)))
  RichPeak2D p;
  p.setIntensity(123.456f);
  p.setMetaValue("cluster_id",4711);
  
  RichPeak2D copy_of_p(p);

  TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456f)
  TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"),DataValue(4711));
END_SECTION

START_SECTION((RichPeak2D(RichPeak2D &&rhs)))
{
  // Ensure that RichPeak2D has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(RichPeak2D(std::declval<RichPeak2D&&>())), true)

  Peak2D::PositionType pos;
  pos[0] = 21.21;
  pos[1] = 22.22;
  RichPeak2D p;
  p.setIntensity(123.456f);
  p.setMetaValue("cluster_id",4711);
  p.setPosition(pos);
  
  RichPeak2D copy_of_p(std::move(p));

  TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456f)
  TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"),DataValue(4711));

  auto i2 = copy_of_p.getIntensity();
  auto pos2 = copy_of_p.getPosition();
  TEST_REAL_SIMILAR(i2, 123.456)

  TEST_REAL_SIMILAR(pos2[0], 21.21)
  TEST_REAL_SIMILAR(pos2[1], 22.22)
}
END_SECTION

START_SECTION((RichPeak2D(const Peak2D &p)))
  Peak2D p;
  p.setIntensity(123.456f);
  
  RichPeak2D copy_of_p(p);

  TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456f)
END_SECTION    
    
START_SECTION((explicit RichPeak2D(const PositionType& pos, const IntensityType in)))
  RichPeak2D p(RichPeak2D::PositionType(21.21, 22.22), 123.456f);
  RichPeak2D copy_of_p(p);
  TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456)
  TEST_REAL_SIMILAR(copy_of_p.getPosition()[0], 21.21)
  TEST_REAL_SIMILAR(copy_of_p.getPosition()[1], 22.22)
END_SECTION

START_SECTION((RichPeak2D& operator=(const RichPeak2D &rhs)))
  RichPeak2D p;
  p.setIntensity(123.456f);
  p.setMetaValue("cluster_id",4711);
  
  RichPeak2D copy_of_p;
  copy_of_p = p;

  TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456f)
  TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"),DataValue(4711));
END_SECTION
    
START_SECTION((RichPeak2D& operator=(const Peak2D &rhs)))
  Peak2D p;
  p.setIntensity(123.456f);
  
  RichPeak2D copy_of_p;
  copy_of_p.setMetaValue("cluster_id",4711);
  copy_of_p = p;

  TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456f)
  TEST_EQUAL(copy_of_p.isMetaEmpty(), true);
END_SECTION
    
START_SECTION((bool operator == (const RichPeak2D& rhs) const))
  RichPeak2D p1, p2;
  TEST_TRUE(p1 == p2)
  
  p1.setIntensity(5.0f);
  TEST_EQUAL(p1==p2, false)
  p2.setIntensity(5.0f);
  TEST_TRUE(p1 == p2)

  p1.setMetaValue("cluster_id",4711);
  TEST_EQUAL(p1==p2, false)
  p1.removeMetaValue("cluster_id");
  TEST_TRUE(p1 == p2)    
END_SECTION

START_SECTION((bool operator != (const RichPeak2D& rhs) const))
  RichPeak2D p1, p2;
  TEST_EQUAL(p1!=p2, false)
  
  p1.setIntensity(5.0f);
  TEST_FALSE(p1 == p2)
  p2.setIntensity(5.0f);
  TEST_EQUAL(p1!=p2, false)

  p1.setMetaValue("cluster_id",4711);
  TEST_FALSE(p1 == p2)
  p1.removeMetaValue("cluster_id");
  TEST_EQUAL(p1!=p2, false)  
END_SECTION

START_SECTION(([EXTRA] meta info with copy constructor))
  RichPeak2D p;
  p.setMetaValue(2,String("bla"));
   RichPeak2D p2(p);
  TEST_EQUAL(p.getMetaValue(2), "bla")
  TEST_EQUAL(p2.getMetaValue(2), "bla")
   p.setMetaValue(2,String("bluff"));
  TEST_EQUAL(p.getMetaValue(2), "bluff")
  TEST_EQUAL(p2.getMetaValue(2), "bla")
END_SECTION

START_SECTION(([EXTRA] meta info with assignment))
  RichPeak2D p;
  p.setMetaValue(2,String("bla"));
  RichPeak2D p2 = p;
  TEST_EQUAL(p.getMetaValue(2), "bla")
  TEST_EQUAL(p2.getMetaValue(2), "bla")
  p.setMetaValue(2,String("bluff"));
  TEST_EQUAL(p.getMetaValue(2), "bluff")
  TEST_EQUAL(p2.getMetaValue(2), "bla")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
