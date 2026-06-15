// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Patrick Boschmann$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/IMAGING/MSImagingRegion.h>
#include <algorithm>

namespace OpenMS
{
MSImagingRegion MSImagingRegion::rectangle(Size id, const String& name, UInt min_x, UInt min_y, UInt max_x, UInt max_y)
{
  if (max_x < min_x || max_y < min_y)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Coordinate maximum is less than minimum",
                                  "min=(" + String(min_x) + "," + String(min_y) + ") max=(" + String(max_x) + "," + String(max_y) + ")");
  }

  MSImagingRegion reg;
  reg.id_ = id;
  reg.name_ = name;
  reg.min_x_ = min_x;
  reg.min_y_ = min_y;
  reg.max_x_ = max_x;
  reg.max_y_ = max_y;
  reg.shape_ = Shape::Rectangle;

  return reg;
}

MSImagingRegion MSImagingRegion::fromMask(Size id, const String& name, UInt origin_x, UInt origin_y, UInt width, UInt height, std::vector<bool> mask)
{
  if (width == 0 || height == 0)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "width or height are 0",
                                  "(w,h): (" + String(width) + "," + String(height) + ")");
  }
  if (mask.size() != size_t(width) * height || std::none_of(mask.begin(), mask.end(), [](bool b) { return b; }))
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "invalid mask",
                                  "(w,h): (" + String(width) + "," + String(height) + ")");
  }
  MSImagingRegion reg;
  reg.id_ = id;
  reg.name_ = name;
  reg.min_x_ = origin_x;
  reg.min_y_ = origin_y;
  reg.max_x_ = origin_x + width - 1;
  reg.max_y_ = origin_y + height - 1;
  reg.mask_ = std::move(mask);
  reg.shape_ = Shape::Mask;

  return reg;
}

const String& MSImagingRegion::getName() const
{ return name_; }


MSImagingRegion::Shape MSImagingRegion::getShape() const
{ return shape_; }

Size MSImagingRegion::getId() const
{ return id_; }

UInt MSImagingRegion::getMinX() const
{ return min_x_; }

UInt MSImagingRegion::getMinY() const
{ return min_y_; }

UInt MSImagingRegion::getMaxX() const
{ return max_x_; }

UInt MSImagingRegion::getMaxY() const
{ return max_y_; }

UInt MSImagingRegion::getBBoxWidth() const
{ return max_x_ - min_x_ + 1; }

UInt MSImagingRegion::getBBoxHeight() const
{ return max_y_ - min_y_ + 1; }

const std::vector<bool>& MSImagingRegion::getMask() const
{ return mask_; }

Size MSImagingRegion::area() const
{
  if (shape_ == Shape::Rectangle) { return static_cast<Size>(getBBoxWidth()) * getBBoxHeight(); }
  return static_cast<Size>(std::count(mask_.begin(), mask_.end(), true));
}

bool MSImagingRegion::contains(UInt x, UInt y) const
{
  if (x < min_x_ || x > max_x_ || y < min_y_ || y > max_y_) { return false; }

  if (shape_ == Shape::Rectangle) { return true; }

  // row major indexing, convert to within bbox coordinates
  const Size idx = (y - min_y_) * getBBoxWidth() + (x - min_x_);
  return mask_[idx];
}
} // namespace OpenMS