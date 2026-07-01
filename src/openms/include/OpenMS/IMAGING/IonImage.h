// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/RangeManager.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Dense W x H grid of ion intensities with a per-pixel mask.

    Storage is row-major: index = y * W + x. Pixels are masked-out by
    default; setIntensity marks them present. The m/z window the image
    was extracted from is stored alongside the data for traceability.

    3D MSI is intentionally not modeled here. Serial-section experiments
    are a stack of independent 2D acquisitions and should be handled as a
    collection of IonImage objects (one per section), matching the OpenMS
    "vector of units" idiom.
  */
  class OPENMS_DLLAPI IonImage final
  {
  public:
    /// @brief Default-construct an empty 0 x 0 image.
    IonImage() = default;

    /**
      @brief Constructs and zero-initializes a @p width x @p height image; all pixels invalid.
      @param[in] width  Number of columns.
      @param[in] height Number of rows.
    */
    IonImage(UInt width, UInt height);

    /**
      @brief Resizes and zero-initializes; all pixels become invalid.
      @param[in] width  Number of columns.
      @param[in] height Number of rows.
    */
    void resize(UInt width, UInt height);

    /// @brief Image width.
    /// @return Number of columns.
    UInt getWidth() const;

    /// @brief Image height.
    /// @return Number of rows.
    UInt getHeight() const;

    /**
      @brief True once setIntensity has been called for the cell at (@p x, @p y).
      @param[in] x Column index.
      @param[in] y Row index.
      @return false for out-of-bounds coordinates or cells never written.
    */
    bool hasPixel(UInt x, UInt y) const;

    /**
      @brief Intensity at (@p x, @p y); 0.0 if never set.
      @param[in] x Column index.
      @param[in] y Row index.
      @return Stored intensity, or 0.0 if hasPixel(@p x, @p y) is false.
      @throws Exception::IndexOverflow on out-of-bounds access.
    */
    double getIntensity(UInt x, UInt y) const;

    /**
      @brief Stores @p intensity at (@p x, @p y) and marks the cell valid.
      @param[in] x         Column index.
      @param[in] y         Row index.
      @param[in] intensity Value to store.
      @throws Exception::IndexOverflow on out-of-bounds access.
    */
    void setIntensity(UInt x, UInt y, double intensity);

    /**
      @brief Records the m/z window the image was extracted from.
      @param[in] range m/z window.
    */
    void setMzRange(const RangeMZ& range);

    /// @brief m/z window the image was extracted from.
    /// @return Stored range; empty if never set.
    const RangeMZ& getMzRange() const;

    /// @brief Raw row-major intensity buffer.
    /// @return Reference to the internal buffer of size width * height.
    const std::vector<double>& getData() const;

    /// @brief Parallel pixel mask (same indexing as getData()).
    /// @return Reference to the mask of size width * height.
    const std::vector<bool>& getMask() const;

  private:
    UInt width_ = 0;
    UInt height_ = 0;
    std::vector<double> intensities_;
    std::vector<bool> mask_;
    RangeMZ mz_range_;

    Size linearIndex_(UInt x, UInt y) const;
  };

} // namespace OpenMS
