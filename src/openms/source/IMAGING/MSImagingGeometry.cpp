// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Patrick Boschmann $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/IMAGING/MSImagingGeometry.h>
#include <algorithm>

namespace OpenMS
{

void MSImagingGeometry::setDimensions(UInt width, UInt height)
{
  width_ = width;
  height_ = height;
}

UInt MSImagingGeometry::getWidth() const
{ return width_; }

UInt MSImagingGeometry::getHeight() const
{ return height_; }

void MSImagingGeometry::setPixelSize(double x, double y, const String& unit)
{
  pixel_size_x_ = x;
  pixel_size_y_ = y;
  pixel_size_unit_ = unit;
}

double MSImagingGeometry::getPixelSizeX() const
{ return pixel_size_x_; }

double MSImagingGeometry::getPixelSizeY() const
{ return pixel_size_y_; }

const String& MSImagingGeometry::getPixelSizeUnit() const
{ return pixel_size_unit_; }

void MSImagingGeometry::addPixel(UInt x, UInt y, Size spectrum_index)
{
  if (width_ > 0 && height_ > 0 && (x >= width_ || y >= height_))
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Pixel coordinate outside configured geometry bounds",
                                  String(x) + "," + String(y));
  }

  const UInt64 key = packKey_(x, y);
  if (lookup_.contains(key))
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Duplicate pixel coordinate", String(x) + "," + String(y));
  }

  lookup_.emplace(key, pixels_.size());
  pixels_.push_back(Pixel {x, y, spectrum_index});
}

bool MSImagingGeometry::hasPixel(UInt x, UInt y) const
{ return lookup_.contains(packKey_(x, y)); }

Size MSImagingGeometry::getSpectrumIndex(UInt x, UInt y) const
{
  auto it = lookup_.find(packKey_(x, y));
  if (it == lookup_.end()) { throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(x) + "," + String(y)); }
  return pixels_[it->second].spectrum_index;
}

const std::vector<MSImagingGeometry::Pixel>& MSImagingGeometry::getPixels() const
{ return pixels_; }

Size MSImagingGeometry::getNumberOfPixels() const
{ return pixels_.size(); }

Size MSImagingGeometry::getNumberOfRegions() const
{ return regions_.size(); }

Size MSImagingGeometry::regionOf(UInt x, UInt y) const
{
  if (! hasPixel(x, y)) return NO_REGION;
  const auto& regions = getRegions();
  for (const auto& region : regions)
  {
    if (region.contains(x, y)) return region.getId();
  }
  return NO_REGION;
}

std::vector<Size> MSImagingGeometry::getRegionPixels(Size id) const
{
  std::vector<Size> result;
  const auto& region = getRegion(id);
  const auto& pixels = getPixels();
  for (Size i = 0; i < pixels.size(); ++i)
  {
    if (region.contains(pixels[i].x, pixels[i].y)) { result.push_back(i); }
  }
  return result;
}

std::vector<Size> MSImagingGeometry::getRegionSpectrumIndices(Size id) const
{
  std::vector<Size> result;
  const auto& region = getRegion(id);
  const auto& pixels = getPixels();
  for (const auto& pixel : pixels)
  {
    if (region.contains(pixel.x, pixel.y)) { result.push_back(pixel.spectrum_index); }
  }
  return result;
}

const std::vector<MSImagingRegion>& MSImagingGeometry::getRegions() const
{ return regions_; }

const MSImagingRegion& MSImagingGeometry::getRegion(Size id) const
{
  auto it = region_id_to_index_.find(id);
  if (it == region_id_to_index_.end()) { throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(id)); }
  return regions_[it->second];
}

void MSImagingGeometry::clearRegions()
{
  regions_.clear();
  region_id_to_index_.clear();
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
  clearRegions();
}

void MSImagingGeometry::addRegion(const MSImagingRegion& region)
{
  if (region.getId() == NO_REGION)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "no region is not a valid region", String(NO_REGION));
  }
  if (region_id_to_index_.contains(region.getId()))
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Duplicate region ID", String(region.getId()));
  }
  const auto& pixels = getPixels();
  for (const auto& pixel : pixels)
  {
    if (! region.contains(pixel.x, pixel.y)) continue;
    const UInt px = pixel.x, py = pixel.y;
    if (std::any_of(regions_.begin(), regions_.end(), [px, py](const MSImagingRegion& r) { return r.contains(px, py); }))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Regions must be disjoint", String(px) + "," + String(py));
    }
  }
  region_id_to_index_[region.getId()] = regions_.size();
  regions_.push_back(region);
}

void MSImagingGeometry::removeRegion(Size id)
{
  auto it = region_id_to_index_.find(id);
  if (it == region_id_to_index_.end()) { throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(id)); }
  const Size region_index = it->second;
  regions_.erase(std::next(regions_.begin(), static_cast<ptrdiff_t>(region_index)));
  region_id_to_index_.clear();
  for (Size i = 0; i < regions_.size(); ++i)
  {
    region_id_to_index_[regions_[i].getId()] = i;
  }
}

UInt64 MSImagingGeometry::packKey_(UInt x, UInt y)
{ return (static_cast<UInt64>(y) << 32) | static_cast<UInt64>(x); }

constexpr Size MSImagingGeometry::NO_REGION; // define const out of line for nanobind compat

} // namespace OpenMS
