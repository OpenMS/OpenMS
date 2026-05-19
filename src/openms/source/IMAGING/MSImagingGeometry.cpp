// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/IMAGING/MSImagingGeometry.h>

#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

  namespace
  {
    constexpr UInt64 MAX_AXIS_VALUE = static_cast<UInt64>(1) << 21; // 2^21 = 2097152
  }

  void MSImagingGeometry::setDimensions(UInt width, UInt height, UInt depth)
  {
    width_ = width;
    height_ = height;
    depth_ = depth;
  }

  UInt MSImagingGeometry::getWidth() const
  {
    return width_;
  }

  UInt MSImagingGeometry::getHeight() const
  {
    return height_;
  }

  UInt MSImagingGeometry::getDepth() const
  {
    return depth_;
  }

  void MSImagingGeometry::setPixelSize(double x, double y, double z, const String& unit)
  {
    pixel_size_x_ = x;
    pixel_size_y_ = y;
    pixel_size_z_ = z;
    pixel_size_unit_ = unit;
  }

  double MSImagingGeometry::getPixelSizeX() const
  {
    return pixel_size_x_;
  }

  double MSImagingGeometry::getPixelSizeY() const
  {
    return pixel_size_y_;
  }

  double MSImagingGeometry::getPixelSizeZ() const
  {
    return pixel_size_z_;
  }

  const String& MSImagingGeometry::getPixelSizeUnit() const
  {
    return pixel_size_unit_;
  }

  void MSImagingGeometry::addPixel(UInt x, UInt y, Size spectrum_index)
  {
    addPixel(x, y, 0, spectrum_index);
  }

  void MSImagingGeometry::addPixel(UInt x, UInt y, UInt z, Size spectrum_index)
  {
    if (static_cast<UInt64>(x) >= MAX_AXIS_VALUE
        || static_cast<UInt64>(y) >= MAX_AXIS_VALUE
        || static_cast<UInt64>(z) >= MAX_AXIS_VALUE)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Pixel coordinate exceeds bit-packing limit (2^21)",
                                    String(x) + "," + String(y) + "," + String(z));
    }

    const UInt64 key = packKey_(x, y, z);
    if (lookup_.find(key) != lookup_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Duplicate pixel coordinate",
                                    String(x) + "," + String(y) + "," + String(z));
    }

    lookup_.emplace(key, pixels_.size());
    pixels_.push_back(Pixel{x, y, z, spectrum_index});
  }

  bool MSImagingGeometry::hasPixel(UInt x, UInt y, UInt z) const
  {
    if (static_cast<UInt64>(x) >= MAX_AXIS_VALUE
        || static_cast<UInt64>(y) >= MAX_AXIS_VALUE
        || static_cast<UInt64>(z) >= MAX_AXIS_VALUE)
    {
      return false;
    }
    return lookup_.find(packKey_(x, y, z)) != lookup_.end();
  }

  Size MSImagingGeometry::getSpectrumIndex(UInt x, UInt y, UInt z) const
  {
    if (static_cast<UInt64>(x) < MAX_AXIS_VALUE
        && static_cast<UInt64>(y) < MAX_AXIS_VALUE
        && static_cast<UInt64>(z) < MAX_AXIS_VALUE)
    {
      auto it = lookup_.find(packKey_(x, y, z));
      if (it != lookup_.end()) return pixels_[it->second].spectrum_index;
    }
    throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     String(x) + "," + String(y) + "," + String(z));
  }

  const std::vector<MSImagingGeometry::Pixel>& MSImagingGeometry::getPixels() const
  {
    return pixels_;
  }

  Size MSImagingGeometry::getNumberOfPixels() const
  {
    return pixels_.size();
  }

  void MSImagingGeometry::clear()
  {
    width_ = 0;
    height_ = 0;
    depth_ = 1;
    pixel_size_x_ = 1.0;
    pixel_size_y_ = 1.0;
    pixel_size_z_ = 1.0;
    pixel_size_unit_ = "micrometer";
    pixels_.clear();
    lookup_.clear();
  }

  UInt64 MSImagingGeometry::packKey_(UInt x, UInt y, UInt z)
  {
    return (static_cast<UInt64>(z) << 42)
           | (static_cast<UInt64>(y) << 21)
           | static_cast<UInt64>(x);
  }

} // namespace OpenMS
