// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  Timo Sachsenberg$
// $Authors: Patrick Boschmann$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief A spatial region within an MSI dataset, in global pixel coordinates.

    A region is a pure-geometry footprint, either an axis-aligned rectangular
    bounding box (@p Shape::Rectangle) or a per-pixel bitmask within a bounding
    box (@p Shape::Mask). It knows nothing about acquired pixels or spectra and
    is reusable as a bare footprint (e.g. a microscopy annotation). Construct
    via the rectangle() / fromMask() factories. Coordinates are zero-based and
    bounding boxes are inclusive on both ends.
  */
class OPENMS_DLLAPI MSImagingRegion final
{
public:
  /// @brief Discriminates how a region's footprint is represented.
  enum class Shape
  {
    Rectangle, ///< grid aligned bounding box (every cell inside the bbox is part of the region)
    Mask       ///< per pixel bitmask inside the bounding box
  };

  /**
    @brief Creates a rectangular region spanning the inclusive bounding box.
    @param[in] id    Region identifier.
    @param[in] name  Human-readable region name.
    @param[in] min_x Leftmost column (inclusive, zero-based).
    @param[in] min_y Top row (inclusive, zero-based).
    @param[in] max_x Rightmost column (inclusive, zero-based).
    @param[in] max_y Bottom row (inclusive, zero-based).
    @return A region with Shape::Rectangle.
    @throws Exception::InvalidValue if @p max_x < @p min_x or @p max_y < @p min_y.
  */
  static MSImagingRegion rectangle(Size id, const std::string& name, UInt min_x, UInt min_y, UInt max_x, UInt max_y);

  /**
    @brief Creates a masked region from a row-major bitmask over a bounding box.

    The global bounding box is [origin_x, origin_x + width - 1] x
    [origin_y, origin_y + height - 1]; @p mask is stored bbox-local.
    @param[in] id       Region identifier.
    @param[in] name     Human-readable region name.
    @param[in] origin_x Leftmost column of the bounding box (zero-based).
    @param[in] origin_y Top row of the bounding box (zero-based).
    @param[in] width    Bounding box width in pixels.
    @param[in] height   Bounding box height in pixels.
    @param[in] mask     Row-major bitmask (true = inside), size == width * height.
    @return A region with Shape::Mask.
    @throws Exception::InvalidValue if @p width or @p height is 0, if
            @p mask.size() != width * height, or if @p mask is all-false.
  */
  static MSImagingRegion fromMask(Size id, const std::string& name, UInt origin_x, UInt origin_y, UInt width, UInt height, std::vector<bool> mask);

  /// @brief Region name.
  /// @return Reference to the stored name.
  const std::string& getName() const;

  /**
    @brief Tests whether the global coordinate (@p x, @p y) is inside the region.
    @param[in] x Column index (zero-based).
    @param[in] y Row index (zero-based).
    @return true if (@p x, @p y) lies within the footprint (bbox for a Rectangle,
            a set bit for a Mask).
  */
  bool contains(UInt x, UInt y) const;

  /**
    @brief Tests whether this region's footprint geometrically overlaps @p other.

    Symmetric and independent of any acquired pixels: returns true iff some
    global coordinate lies inside both footprints.
    @param[in] other The other region.
    @return true if the footprints share at least one coordinate.
  */
  bool intersects(const MSImagingRegion& other) const;

  /// @brief Region shape (Rectangle or Mask).
  /// @return The shape discriminator.
  Shape getShape() const;

  /// @brief Region identifier.
  /// @return The stored id.
  Size getId() const;

  /// @brief Leftmost column of the global bounding box (inclusive).
  /// @return min x.
  UInt getMinX() const;
  /// @brief Top row of the global bounding box (inclusive).
  /// @return min y.
  UInt getMinY() const;
  /// @brief Rightmost column of the global bounding box (inclusive).
  /// @return max x.
  UInt getMaxX() const;
  /// @brief Bottom row of the global bounding box (inclusive).
  /// @return max y.
  UInt getMaxY() const;

  /// @brief Bounding box width in pixels (max_x - min_x + 1).
  /// @return Bounding box width.
  UInt getBBoxWidth() const;
  /// @brief Bounding box height in pixels (max_y - min_y + 1).
  /// @return Bounding box height.
  UInt getBBoxHeight() const;

  /// @brief Bbox-local row-major bitmask; empty for a Rectangle.
  /// @return Reference to the mask (branch on getShape() before use).
  const std::vector<bool>& getMask() const;

  /// @brief Number of pixels inside the region.
  /// @return Bounding-box area for a Rectangle, set-bit count (popcount) for a Mask.
  Size area() const;

private:
  Size id_ {0};
  Shape shape_ = Shape::Rectangle;
  std::string name_;
  std::vector<bool> mask_; // empty if shape_ is set, otherwise rect is bounding box
  UInt min_x_ {0};
  UInt min_y_ {0};
  UInt max_x_ {0};
  UInt max_y_ {0};
};
} // namespace OpenMS
