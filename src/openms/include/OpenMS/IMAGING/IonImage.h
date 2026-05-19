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
    @brief Dense W x H x D grid of ion intensities with a per-pixel validity mask.

    Storage is row-major, z-major: index = (z * H + y) * W + x.

    Pixels are invalid by default; setIntensity marks them valid. The m/z
    window the image was extracted from is stored alongside the data for
    traceability.
  */
  class OPENMS_DLLAPI IonImage final
  {
  public:
    IonImage() = default;
    IonImage(UInt width, UInt height, UInt depth = 1);

    /// Resize and zero-initialize; all pixels become invalid.
    void resize(UInt width, UInt height, UInt depth = 1);

    /// @name Dimensions
    ///@{
    UInt getWidth() const;
    UInt getHeight() const;
    UInt getDepth() const;
    ///@}

    /// @brief True once setIntensity has been called for the given cell.
    bool hasPixel(UInt x, UInt y, UInt z = 0) const;

    /// @brief Intensity at (x,y,z); 0.0 if never set.
    /// @throws Exception::IndexOverflow on out-of-bounds access.
    double getIntensity(UInt x, UInt y, UInt z = 0) const;

    /// @brief Sets the intensity and marks the cell valid (z = 0).
    /// @throws Exception::IndexOverflow on out-of-bounds access.
    void setIntensity(UInt x, UInt y, double intensity);

    /// @brief Sets the intensity and marks the cell valid.
    /// @throws Exception::IndexOverflow on out-of-bounds access.
    void setIntensity(UInt x, UInt y, UInt z, double intensity);

    /// @name m/z window the image was extracted from
    ///@{
    void setMzRange(const RangeMZ& range);
    const RangeMZ& getMzRange() const;
    ///@}

    /// @brief Raw row-major, z-major intensity buffer.
    const std::vector<double>& getData() const;

    /// @brief Parallel validity mask (same indexing as getData()).
    const std::vector<bool>& getValidity() const;

  private:
    UInt width_ = 0;
    UInt height_ = 0;
    UInt depth_ = 1;
    std::vector<double> intensities_;
    std::vector<bool> valid_pixels_;
    RangeMZ mz_range_;

    Size linearIndex_(UInt x, UInt y, UInt z) const;
  };

} // namespace OpenMS
