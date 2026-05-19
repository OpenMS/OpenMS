// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <unordered_map>
#include <vector>

namespace OpenMS
{
  /**
    @brief Pixel grid metadata and (x,y,z) -> spectrum_index lookup for MSI data.

    Coordinates are zero-based. imzML files are one-based; the Phase 2
    loader (ImzMLFile) is responsible for normalizing those to zero-based
    coordinates before populating this geometry.

    Pixel coordinates are bit-packed into a uint64_t key (21 bits per axis),
    so coordinates must be strictly less than 2^21 = 2097152.
  */
  class OPENMS_DLLAPI MSImagingGeometry final
  {
  public:
    /// @brief A pixel in the imaging grid, linked to one spectrum in the experiment.
    struct Pixel
    {
      UInt x = 0;
      UInt y = 0;
      UInt z = 0;
      Size spectrum_index = 0;
    };

    /// @name Image dimensions (number of pixels along each axis)
    ///@{
    void setDimensions(UInt width, UInt height, UInt depth = 1);
    UInt getWidth() const;
    UInt getHeight() const;
    UInt getDepth() const;
    ///@}

    /// @name Physical pixel size
    ///@{
    void setPixelSize(double x, double y, double z = 1.0,
                      const String& unit = "micrometer");
    double getPixelSizeX() const;
    double getPixelSizeY() const;
    double getPixelSizeZ() const;
    const String& getPixelSizeUnit() const;
    ///@}

    /// @brief Adds a pixel at (x, y, z = 0) bound to @p spectrum_index.
    /// @throws Exception::InvalidValue on duplicate coordinates or when any
    ///         coordinate is >= 2^21.
    void addPixel(UInt x, UInt y, Size spectrum_index);

    /// @brief Adds a pixel at (x, y, z) bound to @p spectrum_index.
    /// @throws Exception::InvalidValue on duplicate coordinates or when any
    ///         coordinate is >= 2^21.
    void addPixel(UInt x, UInt y, UInt z, Size spectrum_index);

    /// @brief True if a pixel exists at (x, y, z).
    bool hasPixel(UInt x, UInt y, UInt z = 0) const;

    /// @brief Returns the spectrum index for (x, y, z).
    /// @throws Exception::ElementNotFound if no pixel exists at that coordinate.
    Size getSpectrumIndex(UInt x, UInt y, UInt z = 0) const;

    /// @brief Pixels in insertion order.
    const std::vector<Pixel>& getPixels() const;

    /// @brief Total number of pixels with a bound spectrum.
    Size getNumberOfPixels() const;

    /// @brief Reset all state (dimensions, pixel size, pixels, lookup).
    void clear();

  private:
    UInt width_ = 0;
    UInt height_ = 0;
    UInt depth_ = 1;
    double pixel_size_x_ = 1.0;
    double pixel_size_y_ = 1.0;
    double pixel_size_z_ = 1.0;
    String pixel_size_unit_ = "micrometer";
    std::vector<Pixel> pixels_;
    std::unordered_map<UInt64, Size> lookup_;

    /// 21 bits per axis; coordinates >= 2^21 must be rejected by callers.
    static UInt64 packKey_(UInt x, UInt y, UInt z);
  };

} // namespace OpenMS
