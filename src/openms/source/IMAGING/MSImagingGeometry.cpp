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

  void MSImagingGeometry::setDimensions(UInt width, UInt height)
  {
    width_ = width;
    height_ = height;
  }

  UInt MSImagingGeometry::getWidth() const
  {
    return width_;
  }

  UInt MSImagingGeometry::getHeight() const
  {
    return height_;
  }

  void MSImagingGeometry::setPixelSize(double x, double y, const String& unit)
  {
    pixel_size_x_ = x;
    pixel_size_y_ = y;
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

  const String& MSImagingGeometry::getPixelSizeUnit() const
  {
    return pixel_size_unit_;
  }

  void MSImagingGeometry::addPixel(UInt x, UInt y, Size spectrum_index)
  {
    if (width_ > 0 && height_ > 0 && (x >= width_ || y >= height_))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Pixel coordinate outside configured geometry bounds",
                                    String(x) + "," + String(y));
    }

    const UInt64 key = packKey_(x, y);
    if (lookup_.find(key) != lookup_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Duplicate pixel coordinate",
                                    String(x) + "," + String(y));
    }

    lookup_.emplace(key, pixels_.size());
    pixels_.push_back(Pixel{x, y, spectrum_index});
  }

  bool MSImagingGeometry::hasPixel(UInt x, UInt y) const
  {
    return lookup_.find(packKey_(x, y)) != lookup_.end();
  }

  Size MSImagingGeometry::getSpectrumIndex(UInt x, UInt y) const
  {
    auto it = lookup_.find(packKey_(x, y));
    if (it == lookup_.end())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       String(x) + "," + String(y));
    }
    return pixels_[it->second].spectrum_index;
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
    pixel_size_x_ = 1.0;
    pixel_size_y_ = 1.0;
    pixel_size_unit_ = "micrometer";
    pixels_.clear();
    lookup_.clear();
  }

  UInt64 MSImagingGeometry::packKey_(UInt x, UInt y)
  {
    return (static_cast<UInt64>(y) << 32) | static_cast<UInt64>(x);
  }

} // namespace OpenMS
