// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
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
  // Default is an empty single-cell Rectangle at the origin.
  TEST_EQUAL(ptr->getShape() == MSImagingRegion::Shape::Rectangle, true)
  TEST_EQUAL(ptr->getId(), 0u)
  TEST_EQUAL(ptr->getMinX(), 0u)
  TEST_EQUAL(ptr->getMinY(), 0u)
  TEST_EQUAL(ptr->getMaxX(), 0u)
  TEST_EQUAL(ptr->getMaxY(), 0u)
  TEST_EQUAL(ptr->contains(0, 0), true)
  TEST_EQUAL(ptr->contains(1, 0), false)
}
END_SECTION

START_SECTION((~MSImagingRegion()))
{
  delete ptr;
}
END_SECTION

START_SECTION((static MSImagingRegion rectangle(Size id, const String& name, UInt min_x, UInt min_y, UInt max_x, UInt max_y)))
{
  MSImagingRegion r = MSImagingRegion::rectangle(7, "tissue A", 2, 3, 5, 8);
  TEST_EQUAL(r.getShape() == MSImagingRegion::Shape::Rectangle, true)
  TEST_EQUAL(r.getId(), 7u)
  TEST_EQUAL(r.getName(), "tissue A")
  TEST_EQUAL(r.getMinX(), 2u)
  TEST_EQUAL(r.getMinY(), 3u)
  TEST_EQUAL(r.getMaxX(), 5u)
  TEST_EQUAL(r.getMaxY(), 8u)
  TEST_EQUAL(r.getBBoxWidth(), 4u)
  TEST_EQUAL(r.getBBoxHeight(), 6u)
  TEST_EQUAL(r.area(), 24u)
  // Mask is empty for a Rectangle.
  TEST_EQUAL(r.getMask().empty(), true)

  // contains() boundaries (closed interval).
  TEST_EQUAL(r.contains(2, 3), true)   // top-left corner
  TEST_EQUAL(r.contains(5, 8), true)   // bottom-right corner
  TEST_EQUAL(r.contains(3, 5), true)   // interior
  TEST_EQUAL(r.contains(1, 3), false)  // left of bbox
  TEST_EQUAL(r.contains(6, 8), false)  // right of bbox
  TEST_EQUAL(r.contains(2, 2), false)  // above bbox
  TEST_EQUAL(r.contains(5, 9), false)  // below bbox

  // single-cell rectangle has area 1
  MSImagingRegion single = MSImagingRegion::rectangle(1, "dot", 4, 4, 4, 4);
  TEST_EQUAL(single.area(), 1u)
  TEST_EQUAL(single.getBBoxWidth(), 1u)
  TEST_EQUAL(single.getBBoxHeight(), 1u)
  TEST_EQUAL(single.contains(4, 4), true)
  TEST_EQUAL(single.contains(4, 5), false)

  // min > max rejected in either dimension
  TEST_EXCEPTION(Exception::InvalidValue, MSImagingRegion::rectangle(1, "bad", 5, 0, 2, 0))
  TEST_EXCEPTION(Exception::InvalidValue, MSImagingRegion::rectangle(1, "bad", 0, 5, 0, 2))
}
END_SECTION

START_SECTION((static MSImagingRegion fromMask(Size id, const String& name, UInt origin_x, UInt origin_y, UInt width, UInt height, std::vector<bool> mask)))
{
  // 3x2 mask (width=3, height=2), origin at (10, 20):
  //   row 0: 1 0 1
  //   row 1: 0 1 0
  std::vector<bool> mask = {true, false, true,
                            false, true, false};
  MSImagingRegion r = MSImagingRegion::fromMask(3, "field", 10, 20, 3, 2, mask);
  TEST_EQUAL(r.getShape() == MSImagingRegion::Shape::Mask, true)
  TEST_EQUAL(r.getId(), 3u)
  TEST_EQUAL(r.getName(), "field")
  TEST_EQUAL(r.getMinX(), 10u)
  TEST_EQUAL(r.getMinY(), 20u)
  TEST_EQUAL(r.getMaxX(), 12u)
  TEST_EQUAL(r.getMaxY(), 21u)
  TEST_EQUAL(r.getBBoxWidth(), 3u)
  TEST_EQUAL(r.getBBoxHeight(), 2u)
  // area is the popcount, not the bbox area
  TEST_EQUAL(r.area(), 3u)
  TEST_EQUAL(r.getMask().size(), 6u)

  // contains() against set bits (global coords)
  TEST_EQUAL(r.contains(10, 20), true)  // local (0,0) set
  TEST_EQUAL(r.contains(12, 20), true)  // local (2,0) set
  TEST_EQUAL(r.contains(11, 21), true)  // local (1,1) set
  // contains() against unset bits inside the bbox
  TEST_EQUAL(r.contains(11, 20), false) // local (1,0) unset
  TEST_EQUAL(r.contains(10, 21), false) // local (0,1) unset
  TEST_EQUAL(r.contains(12, 21), false) // local (2,1) unset
  // outside the bbox entirely
  TEST_EQUAL(r.contains(9, 20), false)
  TEST_EQUAL(r.contains(13, 21), false)
  TEST_EQUAL(r.contains(10, 22), false)

  // zero-dim rejected
  TEST_EXCEPTION(Exception::InvalidValue, MSImagingRegion::fromMask(1, "bad", 0, 0, 0, 2, std::vector<bool>{}))
  TEST_EXCEPTION(Exception::InvalidValue, MSImagingRegion::fromMask(1, "bad", 0, 0, 2, 0, std::vector<bool>{}))
  // size mismatch rejected
  TEST_EXCEPTION(Exception::InvalidValue, MSImagingRegion::fromMask(1, "bad", 0, 0, 2, 2, std::vector<bool>{true, false, true}))
  // all-false mask rejected (degenerate)
  TEST_EXCEPTION(Exception::InvalidValue, MSImagingRegion::fromMask(1, "bad", 0, 0, 2, 2, std::vector<bool>(4, false)))
}
END_SECTION

START_SECTION((bool contains(UInt x, UInt y) const))
{
  NOT_TESTABLE // tested in the factory sections above
}
END_SECTION

START_SECTION((Shape getShape() const))
{
  NOT_TESTABLE // tested in the factory sections above
}
END_SECTION

START_SECTION((Size getId() const))
{
  NOT_TESTABLE // tested in the factory sections above
}
END_SECTION

START_SECTION((const String& getName() const))
{
  NOT_TESTABLE // tested in the factory sections above
}
END_SECTION

START_SECTION((UInt getMinX() const))
{
  NOT_TESTABLE // tested in the factory sections above
}
END_SECTION

START_SECTION((UInt getMinY() const))
{
  NOT_TESTABLE // tested in the factory sections above
}
END_SECTION

START_SECTION((UInt getMaxX() const))
{
  NOT_TESTABLE // tested in the factory sections above
}
END_SECTION

START_SECTION((UInt getMaxY() const))
{
  NOT_TESTABLE // tested in the factory sections above
}
END_SECTION

START_SECTION((UInt getBBoxWidth() const))
{
  NOT_TESTABLE // tested in the factory sections above
}
END_SECTION

START_SECTION((UInt getBBoxHeight() const))
{
  NOT_TESTABLE // tested in the factory sections above
}
END_SECTION

START_SECTION((const std::vector<bool>& getMask() const))
{
  NOT_TESTABLE // tested in the factory sections above
}
END_SECTION

START_SECTION((Size area() const))
{
  NOT_TESTABLE // tested in the factory sections above
}
END_SECTION

START_SECTION((UInt toLocalX(UInt x) const))
{
  MSImagingRegion r = MSImagingRegion::rectangle(1, "r", 2, 3, 5, 8);
  TEST_EQUAL(r.toLocalX(2), 0u)  // left edge
  TEST_EQUAL(r.toLocalX(5), 3u)  // right edge
  TEST_EQUAL(r.toLocalX(4), 2u)  // interior
  // out of bbox in x rejected
  TEST_EXCEPTION(Exception::IndexOverflow, r.toLocalX(1))
  TEST_EXCEPTION(Exception::IndexOverflow, r.toLocalX(6))
}
END_SECTION

START_SECTION((UInt toLocalY(UInt y) const))
{
  MSImagingRegion r = MSImagingRegion::rectangle(1, "r", 2, 3, 5, 8);
  TEST_EQUAL(r.toLocalY(3), 0u)  // top edge
  TEST_EQUAL(r.toLocalY(8), 5u)  // bottom edge
  TEST_EQUAL(r.toLocalY(5), 2u)  // interior
  // out of bbox in y rejected
  TEST_EXCEPTION(Exception::IndexOverflow, r.toLocalY(2))
  TEST_EXCEPTION(Exception::IndexOverflow, r.toLocalY(9))
}
END_SECTION

START_SECTION((Size localIndex(UInt x, UInt y) const))
{
  // bbox 4 wide (x in [2,5]) x 6 tall (y in [3,8]); index = ly * 4 + lx
  MSImagingRegion r = MSImagingRegion::rectangle(1, "r", 2, 3, 5, 8);
  TEST_EQUAL(r.localIndex(2, 3), 0u)   // local (0,0)
  TEST_EQUAL(r.localIndex(5, 3), 3u)   // local (3,0)
  TEST_EQUAL(r.localIndex(2, 4), 4u)   // local (0,1)
  TEST_EQUAL(r.localIndex(5, 8), 23u)  // local (3,5) -> 5*4+3
  // index matches getMask() indexing for a Mask region
  std::vector<bool> mask = {true, false, true,
                            false, true, false};
  MSImagingRegion m = MSImagingRegion::fromMask(2, "m", 10, 20, 3, 2, mask);
  TEST_EQUAL(m.localIndex(10, 20), 0u)  // set bit
  TEST_EQUAL(m.localIndex(12, 20), 2u)  // set bit
  TEST_EQUAL(m.localIndex(11, 21), 4u)  // set bit (1*3+1)
  TEST_EQUAL(m.getMask()[m.localIndex(11, 21)], true)
  TEST_EQUAL(m.getMask()[m.localIndex(11, 20)], false)
  // localIndex ignores membership: an unset Mask cell inside the bbox still maps
  TEST_EQUAL(m.localIndex(10, 21), 3u)
  // out of bbox rejected (either dimension)
  TEST_EXCEPTION(Exception::IndexOverflow, m.localIndex(9, 20))
  TEST_EXCEPTION(Exception::IndexOverflow, m.localIndex(10, 22))
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
