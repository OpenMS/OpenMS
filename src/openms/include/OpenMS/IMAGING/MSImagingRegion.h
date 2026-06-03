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

#include <vector>

namespace OpenMS
{
  /**
    @brief A region footprint in global pixel coordinates for MSI data.

    A region describes a contiguous or arbitrary subset of an imaging-MS
    acquisition (e.g. a separate tissue area / ablation field, or a
    user-drawn selection added at analysis time). It is a pure geometric
    value object: it knows nothing about acquired pixels or spectra and can
    be reused as a bare footprint (e.g. a microscopy annotation).

    Two shapes are supported:
      - @c Rectangle: an axis-aligned box [min_x, max_x] x [min_y, max_y]
        (closed interval; both margins included).
      - @c Mask: an arbitrary bitmask. The bounding box is stored in global
        coordinates; the bitmask itself is stored row-major and
        @b bbox-local for compactness (mirrors the IonImage mask idiom).

    Coordinates are zero-based and global throughout the public API
    (contains() and the bounding-box accessors). For per-region processing
    that works on a bbox-cropped buffer (e.g. a region-local IonImage), the
    toLocalX() / toLocalY() / localIndex() helpers translate a global pixel
    into the region's bbox-local frame. Construct via the rectangle() /
    fromMask() factories; the default constructor yields an empty 0..0
    Rectangle.
  */
  class OPENMS_DLLAPI MSImagingRegion final
  {
  public:
    /// @brief Region footprint shape.
    enum class Shape
    {
      Rectangle, ///< Axis-aligned box, no per-cell mask.
      Mask       ///< Arbitrary bitmask stored bbox-local.
    };

    /// @brief Default-construct an empty single-cell Rectangle at the origin (id 0).
    MSImagingRegion() = default;

    /**
      @brief Builds an axis-aligned rectangular region (global coords, closed interval).
      @param[in] id    Region identifier (unique within a geometry).
      @param[in] name  Human-readable label.
      @param[in] min_x Left column (inclusive).
      @param[in] min_y Top row (inclusive).
      @param[in] max_x Right column (inclusive).
      @param[in] max_y Bottom row (inclusive).
      @return The constructed region.
      @throws Exception::InvalidValue if @p min_x > @p max_x or @p min_y > @p max_y.
    */
    static MSImagingRegion rectangle(Size id, const String& name,
                                     UInt min_x, UInt min_y, UInt max_x, UInt max_y);

    /**
      @brief Builds a bitmask region; the mask is row-major and bbox-local.

      The global bounding box is
      [origin_x, origin_x + width - 1] x [origin_y, origin_y + height - 1].
      Mask element (lx, ly) maps to global (origin_x + lx, origin_y + ly)
      and is stored at index ly * width + lx.

      @param[in] id       Region identifier (unique within a geometry).
      @param[in] name     Human-readable label.
      @param[in] origin_x Global x of the mask's top-left corner.
      @param[in] origin_y Global y of the mask's top-left corner.
      @param[in] width    Mask width in cells (> 0).
      @param[in] height   Mask height in cells (> 0).
      @param[in] mask     Row-major bbox-local bitmask of size width * height.
      @return The constructed region.
      @throws Exception::InvalidValue if @p width or @p height is 0, if
              @p mask.size() != width * height, or if @p mask is all-false
              (degenerate / empty region).
    */
    static MSImagingRegion fromMask(Size id, const String& name,
                                    UInt origin_x, UInt origin_y,
                                    UInt width, UInt height, std::vector<bool> mask);

    /**
      @brief Tests whether the global pixel (@p x, @p y) belongs to this region.
      @param[in] x Global column index.
      @param[in] y Global row index.
      @return true if (@p x, @p y) is inside the bounding box and, for a Mask
              region, the corresponding bit is set.
    */
    bool contains(UInt x, UInt y) const;

    /// @brief Region shape.
    /// @return Rectangle or Mask.
    Shape getShape() const;

    /// @brief Region identifier.
    /// @return Stored id.
    Size getId() const;

    /// @brief Region label.
    /// @return Reference to the stored name.
    const String& getName() const;

    /// @brief Left column of the global bounding box (inclusive).
    UInt getMinX() const;
    /// @brief Top row of the global bounding box (inclusive).
    UInt getMinY() const;
    /// @brief Right column of the global bounding box (inclusive).
    UInt getMaxX() const;
    /// @brief Bottom row of the global bounding box (inclusive).
    UInt getMaxY() const;

    /// @brief Bounding-box width.
    /// @return max_x - min_x + 1.
    UInt getBBoxWidth() const;

    /// @brief Bounding-box height.
    /// @return max_y - min_y + 1.
    UInt getBBoxHeight() const;

    /**
      @brief Raw row-major bbox-local bitmask.
      @return Reference to the mask of size bbox_width * bbox_height for a
              Mask region; an empty vector for a Rectangle region (callers
              use getShape() / contains() instead).
    */
    const std::vector<bool>& getMask() const;

    /**
      @brief Number of cells the region covers.
      @return Rectangle: bbox_width * bbox_height. Mask: number of set bits.
    */
    Size area() const;

    /**
      @brief Translates a global column to the region's bbox-local frame.
      @param[in] x Global column index (must lie within the bounding box).
      @return @p x - min_x (local column, zero-based at the bbox left edge).
      @throws Exception::IndexOverflow if @p x is outside [min_x, max_x].
    */
    UInt toLocalX(UInt x) const;

    /**
      @brief Translates a global row to the region's bbox-local frame.
      @param[in] y Global row index (must lie within the bounding box).
      @return @p y - min_y (local row, zero-based at the bbox top edge).
      @throws Exception::IndexOverflow if @p y is outside [min_y, max_y].
    */
    UInt toLocalY(UInt y) const;

    /**
      @brief Row-major bbox-local index of the global pixel (@p x, @p y).

      Matches the indexing used by getMask() and by a bbox-cropped IonImage:
      index = toLocalY(y) * bbox_width + toLocalX(x). Membership in a Mask
      region is @e not consulted; only the bounding box is checked.

      @param[in] x Global column index (must lie within the bounding box).
      @param[in] y Global row index (must lie within the bounding box).
      @return Local linear index in [0, bbox_width * bbox_height).
      @throws Exception::IndexOverflow if (@p x, @p y) is outside the bounding box.
    */
    Size localIndex(UInt x, UInt y) const;

  private:
    Shape shape_ = Shape::Rectangle;
    Size id_ = 0;
    String name_;
    UInt min_x_ = 0;
    UInt min_y_ = 0;
    UInt max_x_ = 0;
    UInt max_y_ = 0;
    std::vector<bool> mask_; ///< Row-major, bbox-local; empty for Rectangle.
  };

} // namespace OpenMS
