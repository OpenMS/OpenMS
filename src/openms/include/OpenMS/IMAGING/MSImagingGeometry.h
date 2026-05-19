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
    @brief Pixel grid metadata and (x, y) -> spectrum_index lookup for MSI data.

    Coordinates are zero-based. imzML files are one-based; the Phase 2 loader
    (ImzMLFile) is responsible for normalizing those to zero-based coordinates
    before populating this geometry.

    3D MSI is intentionally not modeled here. Serial-section experiments
    should be handled as a collection of MSImagingExperiment objects (one
    per section).
  */
  class OPENMS_DLLAPI MSImagingGeometry final
  {
  public:
    /// @brief A pixel in the imaging grid, linked to one spectrum in the experiment.
    struct Pixel
    {
      UInt x = 0;              ///< Column index (zero-based).
      UInt y = 0;              ///< Row index (zero-based).
      Size spectrum_index = 0; ///< Index into the bound MSExperiment.
    };

    /**
      @brief Sets the image dimensions.
      @param[in] width  Number of columns.
      @param[in] height Number of rows.
    */
    void setDimensions(UInt width, UInt height);

    /// @brief Image width.
    /// @return Number of columns.
    UInt getWidth() const;

    /// @brief Image height.
    /// @return Number of rows.
    UInt getHeight() const;

    /**
      @brief Records the physical pixel size and its unit.
      @param[in] x    Pixel size along x.
      @param[in] y    Pixel size along y.
      @param[in] unit Length unit (default "micrometer").
    */
    void setPixelSize(double x, double y, const String& unit = "micrometer");

    /// @brief Physical pixel size along x.
    /// @return Stored x pixel size.
    double getPixelSizeX() const;

    /// @brief Physical pixel size along y.
    /// @return Stored y pixel size.
    double getPixelSizeY() const;

    /// @brief Unit for the pixel size.
    /// @return Reference to the unit string.
    const String& getPixelSizeUnit() const;

    /**
      @brief Adds a pixel at (@p x, @p y) bound to @p spectrum_index.
      @param[in] x              Column index (zero-based).
      @param[in] y              Row index (zero-based).
      @param[in] spectrum_index Index into the bound MSExperiment.
      @throws Exception::InvalidValue on duplicate coordinates, or if
              dimensions have been set and (@p x, @p y) is outside
              [0, width) x [0, height).
    */
    void addPixel(UInt x, UInt y, Size spectrum_index);

    /**
      @brief Tests pixel presence at (@p x, @p y).
      @param[in] x Column index.
      @param[in] y Row index.
      @return true if a pixel was inserted at that coordinate.
    */
    bool hasPixel(UInt x, UInt y) const;

    /**
      @brief Looks up the spectrum index at (@p x, @p y).
      @param[in] x Column index.
      @param[in] y Row index.
      @return spectrum_index recorded for that pixel.
      @throws Exception::ElementNotFound if no pixel exists at that coordinate.
    */
    Size getSpectrumIndex(UInt x, UInt y) const;

    /// @brief Pixels in insertion order.
    /// @return Reference to the internal pixel list.
    const std::vector<Pixel>& getPixels() const;

    /// @brief Total number of pixels with a bound spectrum.
    /// @return Size of the pixel list.
    Size getNumberOfPixels() const;

    /// @brief Resets all state (dimensions, pixel size, pixels, lookup).
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
