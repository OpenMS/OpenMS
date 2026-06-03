// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/IMAGING/MSImagingRegion.h>

#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

  MSImagingRegion MSImagingRegion::rectangle(Size id, const String& name,
                                             UInt min_x, UInt min_y, UInt max_x, UInt max_y)
  {
    if (min_x > max_x || min_y > max_y)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Rectangle region requires min <= max in both dimensions",
                                    "min=(" + String(min_x) + "," + String(min_y) + "), max=("
                                    + String(max_x) + "," + String(max_y) + ")");
    }

    MSImagingRegion r;
    r.shape_ = Shape::Rectangle;
    r.id_ = id;
    r.name_ = name;
    r.min_x_ = min_x;
    r.min_y_ = min_y;
    r.max_x_ = max_x;
    r.max_y_ = max_y;
    // mask_ stays empty for Rectangle.
    return r;
  }

  MSImagingRegion MSImagingRegion::fromMask(Size id, const String& name,
                                            UInt origin_x, UInt origin_y,
                                            UInt width, UInt height, std::vector<bool> mask)
  {
    if (width == 0 || height == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Mask region requires width > 0 and height > 0",
                                    String(width) + "x" + String(height));
    }
    if (mask.size() != static_cast<size_t>(width) * static_cast<size_t>(height))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Mask size does not match width * height",
                                    String(mask.size()) + " != " + String(width) + "*" + String(height));
    }
    // Reject a degenerate (empty) region: at least one bit must be set.
    bool any_set = false;
    for (const bool b : mask)
    {
      if (b) { any_set = true; break; }
    }
    if (!any_set)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Mask region must have at least one set cell", name);
    }

    MSImagingRegion r;
    r.shape_ = Shape::Mask;
    r.id_ = id;
    r.name_ = name;
    r.min_x_ = origin_x;
    r.min_y_ = origin_y;
    r.max_x_ = origin_x + width - 1;
    r.max_y_ = origin_y + height - 1;
    r.mask_ = std::move(mask);
    return r;
  }

  bool MSImagingRegion::contains(UInt x, UInt y) const
  {
    if (x < min_x_ || x > max_x_ || y < min_y_ || y > max_y_)
    {
      return false;
    }
    if (shape_ == Shape::Rectangle)
    {
      return true;
    }
    // Mask: translate to bbox-local coordinates and test the bit.
    const Size lx = static_cast<Size>(x) - static_cast<Size>(min_x_);
    const Size ly = static_cast<Size>(y) - static_cast<Size>(min_y_);
    const Size idx = ly * static_cast<Size>(getBBoxWidth()) + lx;
    return mask_[idx];
  }

  MSImagingRegion::Shape MSImagingRegion::getShape() const
  {
    return shape_;
  }

  Size MSImagingRegion::getId() const
  {
    return id_;
  }

  const String& MSImagingRegion::getName() const
  {
    return name_;
  }

  UInt MSImagingRegion::getMinX() const
  {
    return min_x_;
  }

  UInt MSImagingRegion::getMinY() const
  {
    return min_y_;
  }

  UInt MSImagingRegion::getMaxX() const
  {
    return max_x_;
  }

  UInt MSImagingRegion::getMaxY() const
  {
    return max_y_;
  }

  UInt MSImagingRegion::getBBoxWidth() const
  {
    return max_x_ - min_x_ + 1;
  }

  UInt MSImagingRegion::getBBoxHeight() const
  {
    return max_y_ - min_y_ + 1;
  }

  const std::vector<bool>& MSImagingRegion::getMask() const
  {
    return mask_;
  }

  Size MSImagingRegion::area() const
  {
    if (shape_ == Shape::Rectangle)
    {
      return static_cast<Size>(getBBoxWidth()) * static_cast<Size>(getBBoxHeight());
    }
    // Mask: count set bits.
    Size count = 0;
    for (const bool b : mask_)
    {
      if (b) { ++count; }
    }
    return count;
  }

} // namespace OpenMS
