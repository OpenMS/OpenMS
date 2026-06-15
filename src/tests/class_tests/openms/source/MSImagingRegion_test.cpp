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

START_SECTION(area() & contains()) {
  NOT_TESTABLE // covered by all region variations
} END_SECTION

  END_TEST
