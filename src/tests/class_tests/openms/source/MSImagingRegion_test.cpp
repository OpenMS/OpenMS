// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Patrick Boschmann $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/IMAGING/MSImagingRegion.h>
///////////////////////////


#include <OpenMS/CONCEPT/Exception.h>


using namespace OpenMS;
using namespace std;

START_TEST(MSImagingRegion, "$Id$")

/////////////////////////////////////////////////////////////

MSImagingRegion* ptr = nullptr;
MSImagingRegion* null_ptr = nullptr;

START_SECTION((MSImagingRegion()))
{
  ptr = new MSImagingRegion();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((~MSImagingRegion()))
{ delete ptr; }
END_SECTION

START_SECTION((MSImagingRegion() defaults))
{
  MSImagingRegion r {};
  TEST_EQUAL(r.getId(), 0u)
  TEST_EQUAL(r.getBBoxWidth(), 1u)
  TEST_EQUAL(r.getBBoxHeight(), 1u)
  TEST_EQUAL(r.getMinX(), 0u)
  TEST_EQUAL(r.getMinY(), 0u)
  TEST_EQUAL(r.getMaxX(), 0u)
  TEST_EQUAL(r.getMaxY(), 0u)
  TEST_EQUAL(r.area(), 1u)
  TEST_EQUAL(r.getMask().size(), 0u)
  TEST_EQUAL(r.getShape(), MSImagingRegion::Shape::Rectangle)
}
END_SECTION

START_SECTION((MSImagingRegion() rectangle))
{
  auto rectangle = MSImagingRegion::rectangle(7, "rect", 2, 3, 5, 6);
  TEST_EQUAL(rectangle.getId(), 7u)
  TEST_EQUAL(rectangle.getName(), "rect")
  TEST_EQUAL(rectangle.getBBoxWidth(), 4u)
  TEST_EQUAL(rectangle.getBBoxHeight(), 4u)
  TEST_EQUAL(rectangle.getMinX(), 2u)
  TEST_EQUAL(rectangle.getMinY(), 3u)
  TEST_EQUAL(rectangle.getMaxX(), 5u)
  TEST_EQUAL(rectangle.getMaxY(), 6u)
  TEST_EQUAL(rectangle.area(), 16u)
  TEST_EQUAL(rectangle.getMask().size(), 0u)
  TEST_EQUAL(rectangle.getShape(), MSImagingRegion::Shape::Rectangle)
  TEST_TRUE(rectangle.contains(2, 3))
  TEST_FALSE(rectangle.contains(1, 4))

  auto single_cell = MSImagingRegion::rectangle(8, "single_cell", 1, 1, 1, 1);
  TEST_EQUAL(single_cell.area(), 1u)
  TEST_TRUE(single_cell.contains(1, 1))
  TEST_FALSE(single_cell.contains(1, 2))
  TEST_FALSE(single_cell.contains(2, 1))

  TEST_EXCEPTION(Exception::InvalidValue, MSImagingRegion::rectangle(9, "exc", 10, 1, 1, 1));
}
END_SECTION

START_SECTION((MSImagingRegion() fromMask))
{
  std::vector<bool> mask {true, false, false, true};
  auto fromMask = MSImagingRegion::fromMask(0, "mask", 10, 20, 2, 2, mask);
  TEST_TRUE(fromMask.getMask() == mask)
  TEST_TRUE(fromMask.getShape() == MSImagingRegion::Shape::Mask)
  TEST_FALSE(fromMask.contains(11, 20))
  TEST_TRUE(fromMask.contains(11, 21))
  TEST_FALSE(fromMask.contains(1, 2))
  TEST_EXCEPTION(Exception::InvalidValue, MSImagingRegion::fromMask(0, "mask", 10, 20, 0, 2, mask));
  TEST_EXCEPTION(Exception::InvalidValue, MSImagingRegion::fromMask(0, "mask", 10, 20, 2, 0, mask));
  TEST_EXCEPTION(Exception::InvalidValue, MSImagingRegion::fromMask(0, "mask", 10, 20, 2, 2,
                                                                    std::vector<bool> {
                                                                      false,
                                                                      false,
                                                                      false,
                                                                      false,
                                                                    }));
  TEST_EXCEPTION(Exception::InvalidValue, MSImagingRegion::fromMask(0, "mask", 10, 20, 2, 2, std::vector<bool> {false, false, false, false, true}));
}
END_SECTION

START_SECTION((bool intersects(const MSImagingRegion& other) const))
{
  //build two overlapping rectangles
  auto a = MSImagingRegion::rectangle(1, "a", 1, 1, 1, 1);
  auto b = MSImagingRegion::rectangle(2, "b", 1, 1, 5, 5);
  // intersects is supposed to be symmetric
  TEST_TRUE(a.intersects(b));
  TEST_TRUE(b.intersects(a));

  auto c = MSImagingRegion::rectangle(3, "c", 10,10, 11,11);
  TEST_FALSE(c.intersects(a));

  //mask with (1,1) and (2,2) set true, overlaps at (1,1)
  auto m1 = MSImagingRegion::fromMask(4, "mask", 1,1,2,2,std::vector<bool>{true, false, false, true});
  TEST_TRUE(m1.intersects(a));

  //overlapping mask with no common pixels (1,2) and (2,1) true
  auto m2 = MSImagingRegion::fromMask(5, "mask2", 1,1,2,2,std::vector<bool>{false, true, true, false});
  TEST_FALSE(m2.intersects(m1));

  //masks with common set:
  auto m3 = MSImagingRegion::fromMask(6, "mask3", 1,1,2,2,std::vector<bool>{true, true, true, false});
  TEST_TRUE(m3.intersects(m2));
}
END_SECTION

START_SECTION(area() & contains()) {
  NOT_TESTABLE // covered by all region variations
} END_SECTION

  END_TEST
