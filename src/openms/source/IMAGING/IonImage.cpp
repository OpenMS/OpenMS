// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/IMAGING/IonImage.h>

#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

  IonImage::IonImage(UInt width, UInt height, UInt depth)
  {
    resize(width, height, depth);
  }

  void IonImage::resize(UInt width, UInt height, UInt depth)
  {
    width_ = width;
    height_ = height;
    depth_ = depth;
    const Size n = static_cast<Size>(width) * static_cast<Size>(height) * static_cast<Size>(depth);
    intensities_.assign(n, 0.0);
    valid_pixels_.assign(n, false);
  }

  UInt IonImage::getWidth() const
  {
    return width_;
  }

  UInt IonImage::getHeight() const
  {
    return height_;
  }

  UInt IonImage::getDepth() const
  {
    return depth_;
  }

  bool IonImage::hasPixel(UInt x, UInt y, UInt z) const
  {
    if (x >= width_ || y >= height_ || z >= depth_) return false;
    return valid_pixels_[linearIndex_(x, y, z)];
  }

  double IonImage::getIntensity(UInt x, UInt y, UInt z) const
  {
    return intensities_[linearIndex_(x, y, z)];
  }

  void IonImage::setIntensity(UInt x, UInt y, double intensity)
  {
    setIntensity(x, y, 0, intensity);
  }

  void IonImage::setIntensity(UInt x, UInt y, UInt z, double intensity)
  {
    const Size idx = linearIndex_(x, y, z);
    intensities_[idx] = intensity;
    valid_pixels_[idx] = true;
  }

  void IonImage::setMzRange(const RangeMZ& range)
  {
    mz_range_ = range;
  }

  const RangeMZ& IonImage::getMzRange() const
  {
    return mz_range_;
  }

  const std::vector<double>& IonImage::getData() const
  {
    return intensities_;
  }

  const std::vector<bool>& IonImage::getValidity() const
  {
    return valid_pixels_;
  }

  Size IonImage::linearIndex_(UInt x, UInt y, UInt z) const
  {
    if (x >= width_ || y >= height_ || z >= depth_)
    {
      const Size attempted = (static_cast<Size>(z) * static_cast<Size>(height_ ? height_ : 1)
                              + static_cast<Size>(y)) * static_cast<Size>(width_ ? width_ : 1)
                             + static_cast<Size>(x);
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     static_cast<SignedSize>(attempted),
                                     intensities_.size());
    }
    return (static_cast<Size>(z) * static_cast<Size>(height_) + static_cast<Size>(y)) * static_cast<Size>(width_) + static_cast<Size>(x);
  }

} // namespace OpenMS
