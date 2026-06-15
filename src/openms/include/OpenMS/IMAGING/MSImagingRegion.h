// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  Timo Sachsenberg$
// $Authors: Patrick Boschmann$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <vector>

namespace OpenMS
{
class OPENMS_DLLAPI MSImagingRegion final
{
public:
  enum class Shape
  {
    Rectangle, // grid aligned bounding box
    Mask       // per pixel mask inside bounding box
  };
  static MSImagingRegion rectangle(Size id, const String& name, UInt min_x, UInt min_y, UInt max_x, UInt max_y);

  static MSImagingRegion fromMask(Size id, const String& name, UInt origin_x, UInt origin_y, UInt width, UInt height, std::vector<bool> mask);

  const String& getName() const;
  bool contains(UInt x, UInt y) const;
  Shape getShape() const;
  Size getId() const;
  UInt getMinX() const;
  UInt getMinY() const;
  UInt getMaxX() const;
  UInt getMaxY() const;
  UInt getBBoxWidth() const;
  UInt getBBoxHeight() const;
  const std::vector<bool>& getMask() const;
  Size area() const; // w*h on a rect, popcount on a mask

private:
  Size id_ {0};
  Shape shape_ = Shape::Rectangle;
  String name_;
  std::vector<bool> mask_; // empty if shape_ is set, otherwise rect is bounding box
  UInt min_x_ {0};
  UInt min_y_ {0};
  UInt max_x_ {0};
  UInt max_y_ {0};
};
} // namespace OpenMS
