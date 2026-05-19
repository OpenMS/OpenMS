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
    @brief Pixel grid metadata and (x,y) -> spectrum_index lookup for MSI data.

    Coordinates are zero-based. imzML files are one-based; the Phase 2
    loader (ImzMLFile) is responsible for normalizing those to zero-based
    coordinates before populating this geometry.

    Pixel coordinates are bit-packed into a uint64_t key (32 bits per
    axis), so any UInt value is admissible.

    3D MSI is intentionally not modeled here. Serial-section experiments
    should be handled as a collection of MSImagingExperiment objects
    (one per section).
  */
  class OPENMS_DLLAPI MSImagingGeometry final
  {
  public:
    /// @brief A pixel in the imaging grid, linked to one spectrum in the experiment.
    struct Pixel
    {
      UInt x = 0;
      UInt y = 0;
      Size spectrum_index = 0;
    };

    /// @name Image dimensions (number of pixels along each axis)
    ///@{
    void setDimensions(UInt width, UInt height);
    UInt getWidth() const;
    UInt getHeight() const;
    ///@}

    /// @name Physical pixel size
    ///@{
    void setPixelSize(double x, double y, const String& unit = "micrometer");
    double getPixelSizeX() const;
    double getPixelSizeY() const;
    const String& getPixelSizeUnit() const;
    ///@}

    /// @brief Adds a pixel at (x, y) bound to @p spectrum_index.
    /// @throws Exception::InvalidValue on duplicate coordinates, or if
    ///         dimensions have been set and (x, y) is outside [0, width) x
    ///         [0, height).
    void addPixel(UInt x, UInt y, Size spectrum_index);

    /// @brief True if a pixel exists at (x, y).
    bool hasPixel(UInt x, UInt y) const;

    /// @brief Returns the spectrum index for (x, y).
    /// @throws Exception::ElementNotFound if no pixel exists at that coordinate.
    Size getSpectrumIndex(UInt x, UInt y) const;

    /// @brief Pixels in insertion order.
    const std::vector<Pixel>& getPixels() const;

    /// @brief Total number of pixels with a bound spectrum.
    Size getNumberOfPixels() const;

    /// @brief Reset all state (dimensions, pixel size, pixels, lookup).
    void clear();

  private:
    UInt width_ = 0;
    UInt height_ = 0;
    double pixel_size_x_ = 1.0;
    double pixel_size_y_ = 1.0;
    String pixel_size_unit_ = "micrometer";
    std::vector<Pixel> pixels_;
    std::unordered_map<UInt64, Size> lookup_;

    static UInt64 packKey_(UInt x, UInt y);
  };

} // namespace OpenMS
