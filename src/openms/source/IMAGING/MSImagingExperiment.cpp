// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Patrick Boschmann $
// --------------------------------------------------------------------------

#include <OpenMS/IMAGING/MSImagingExperiment.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/IMAGING/IonImageExtraction.h>
#include <OpenMS/METADATA/MetaInfoRegistry.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace OpenMS
{

namespace
{
  // Meta value keys resolved once to registry indices: rebuildGeometry() and validate()
  // touch every spectrum, and the by-name accessors would redo the string lookup each time.
  struct PixelMetaIndices
  {
    UInt x;
    UInt y;
    UInt z;
  };

  PixelMetaIndices pixelMetaIndices_()
  {
    // registerName() (not getIndex(), which is UINT_MAX for a name nobody has used yet in this
    // process) returns the existing index when the key is already registered.
    MetaInfoRegistry& reg = MetaInfoInterface::metaRegistry();
    return {reg.registerName(MSImagingExperiment::META_PIXEL_X),
            reg.registerName(MSImagingExperiment::META_PIXEL_Y),
            reg.registerName(MSImagingExperiment::META_PIXEL_Z)};
  }

  // Same contract as MSImagingExperiment::getPixelCoordinate(), on pre-resolved indices.
  bool readPixelCoordinate_(const MSSpectrum& spec, const PixelMetaIndices& idx, UInt& x, UInt& y)
  {
    if (!spec.metaValueExists(idx.x) || !spec.metaValueExists(idx.y))
    {
      return false;
    }
    const Int x_imz = spec.getMetaValue(idx.x);
    const Int y_imz = spec.getMetaValue(idx.y);
    if (x_imz < 1 || y_imz < 1)
    {
      return false;
    }
    if (spec.metaValueExists(idx.z) && static_cast<Int>(spec.getMetaValue(idx.z)) != 1)
    {
      return false;
    }
    x = static_cast<UInt>(x_imz - 1);
    y = static_cast<UInt>(y_imz - 1);
    return true;
  }

  void writePixelCoordinate_(MSSpectrum& spec, const PixelMetaIndices& idx, UInt x, UInt y)
  {
    spec.setMetaValue(idx.x, static_cast<Int>(x) + 1); // 0-based geometry -> 1-based imzML
    spec.setMetaValue(idx.y, static_cast<Int>(y) + 1);
    if (!spec.metaValueExists(idx.z))
    {
      spec.setMetaValue(idx.z, 1);
    }
  }

  std::string coordToString_(UInt x, UInt y)
  {
    return StringUtils::toStr(x) + "," + StringUtils::toStr(y);
  }
} // namespace

MSImagingExperiment::MSImagingExperiment(MSExperiment exp): experiment_(std::move(exp))
{
}

MSImagingExperiment& MSImagingExperiment::operator=(MSExperiment exp)
{
  experiment_ = std::move(exp);
  geometry_.clear();
  return *this;
}

MSExperiment& MSImagingExperiment::getMSExperiment()
{ return experiment_; }

const MSExperiment& MSImagingExperiment::getMSExperiment() const
{ return experiment_; }

void MSImagingExperiment::setMSExperiment(MSExperiment exp)
{ experiment_ = std::move(exp); }

MSImagingGeometry& MSImagingExperiment::getGeometry()
{ return geometry_; }

const MSImagingGeometry& MSImagingExperiment::getGeometry() const
{ return geometry_; }

void MSImagingExperiment::setGeometry(MSImagingGeometry geom)
{
  geometry_ = std::move(geom);
  stampPixelCoordinates_();
}

void MSImagingExperiment::stampPixelCoordinates_()
{
  const PixelMetaIndices idx = pixelMetaIndices_();
  const Size n_spectra = experiment_.getNrSpectra();
  for (const auto& p : geometry_.getPixels())
  {
    if (p.spectrum_index < n_spectra) // dangling pixels are left for validate() to report
    {
      writePixelCoordinate_(experiment_[p.spectrum_index], idx, p.x, p.y);
    }
  }
}

Size MSImagingExperiment::getNumberOfPixels() const
{ return geometry_.getNumberOfPixels(); }

Size MSImagingExperiment::getNrSpectra() const
{ return experiment_.getNrSpectra(); }

Size MSImagingExperiment::size() const
{ return experiment_.getNrSpectra(); }

bool MSImagingExperiment::empty() const
{ return experiment_.empty(); }

MSSpectrum& MSImagingExperiment::getSpectrum(Size index)
{
  if (index >= experiment_.getNrSpectra()) // MSExperiment::getSpectrum() is unchecked
  {
    throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                   static_cast<SignedSize>(index),
                                   static_cast<SignedSize>(experiment_.getNrSpectra()));
  }
  return experiment_[index];
}

const MSSpectrum& MSImagingExperiment::getSpectrum(Size index) const
{
  if (index >= experiment_.getNrSpectra())
  {
    throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                   static_cast<SignedSize>(index),
                                   static_cast<SignedSize>(experiment_.getNrSpectra()));
  }
  return experiment_[index];
}

MSSpectrum& MSImagingExperiment::operator[](Size index)
{ return getSpectrum(index); }

const MSSpectrum& MSImagingExperiment::operator[](Size index) const
{ return getSpectrum(index); }

MSImagingExperiment::Iterator MSImagingExperiment::begin()
{ return experiment_.begin(); }

MSImagingExperiment::Iterator MSImagingExperiment::end()
{ return experiment_.end(); }

MSImagingExperiment::ConstIterator MSImagingExperiment::begin() const
{ return experiment_.begin(); }

MSImagingExperiment::ConstIterator MSImagingExperiment::end() const
{ return experiment_.end(); }

bool MSImagingExperiment::hasPixel(UInt x, UInt y) const
{ return geometry_.hasPixel(x, y); }

void MSImagingExperiment::bindPixel(UInt x, UInt y, Size spectrum_index)
{
  if (spectrum_index >= experiment_.getNrSpectra())
  {
    throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                   static_cast<SignedSize>(spectrum_index),
                                   static_cast<SignedSize>(experiment_.getNrSpectra()));
  }
  geometry_.addPixel(x, y, spectrum_index); // throws on duplicate / out-of-grid before anything is written
  setPixelCoordinate(experiment_[spectrum_index], x, y);
}

  MSSpectrum& MSImagingExperiment::getSpectrum(UInt x, UInt y)
  {
    const Size idx = geometry_.getSpectrumIndex(x, y);
    if (idx >= experiment_.getNrSpectra())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Pixel references missing spectrum",
                                    StringUtils::toStr(idx));
    }
    return experiment_[idx];
  }

  const MSSpectrum& MSImagingExperiment::getSpectrum(UInt x, UInt y) const
  {
    const Size idx = geometry_.getSpectrumIndex(x, y);
    if (idx >= experiment_.getNrSpectra())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Pixel references missing spectrum",
                                    StringUtils::toStr(idx));
    }
    return experiment_[idx];
  }

  IonImage MSImagingExperiment::extractIonImage(double mz, double tolerance_ppm) const
  {
    std::vector<Size> all(geometry_.getNumberOfPixels());
    std::iota(all.begin(), all.end(), Size(0));
    return Internal::extractIonImage(geometry_, mz, tolerance_ppm, all, experiment_.getNrSpectra(),
                                     [this](Size i) -> const MSSpectrum& { return experiment_[i]; });
  }

  IonImage MSImagingExperiment::extractIonImage(double mz, double tolerance_ppm, Size region_id) const
  {
    const auto& indices = geometry_.getRegionPixels(region_id);
    return Internal::extractIonImage(geometry_, mz, tolerance_ppm, indices, experiment_.getNrSpectra(),
                                     [this](Size i) -> const MSSpectrum& { return experiment_[i]; });
  }

  std::vector<Size> MSImagingExperiment::getRegionSpectrumIndices(Size region_id) const
  {
    return geometry_.getRegionSpectrumIndices(region_id);
  }

  bool MSImagingExperiment::getPixelCoordinate(const MSSpectrum& spec, UInt& x, UInt& y)
  {
    return readPixelCoordinate_(spec, pixelMetaIndices_(), x, y);
  }

  void MSImagingExperiment::setPixelCoordinate(MSSpectrum& spec, UInt x, UInt y)
  {
    writePixelCoordinate_(spec, pixelMetaIndices_(), x, y);
  }

  void MSImagingExperiment::bindPixelsFromSpectra(const MSExperiment& exp, MSImagingGeometry& geom)
  {
    const PixelMetaIndices idx = pixelMetaIndices_();
    const UInt width = geom.getWidth();
    const UInt height = geom.getHeight();

    UInt max_x = 0;
    UInt max_y = 0;
    // A broken converter can put every spectrum on the same pixel, so cap what we retain.
    const Size max_listed = 20;
    Size duplicate_count = 0;
    std::vector<Size> duplicate_spectra; // first few spectrum indices dropped from the geometry
    duplicate_spectra.reserve(max_listed);
    const Size n_spectra = exp.getNrSpectra();
    for (Size i = 0; i < n_spectra; ++i)
    {
      const MSSpectrum& spec = exp[i];
      if (!spec.metaValueExists(idx.x) || !spec.metaValueExists(idx.y))
      {
        continue;
      }
      const Int x_imz = spec.getMetaValue(idx.x);
      const Int y_imz = spec.getMetaValue(idx.y);
      if (x_imz < 1 || y_imz < 1)
      {
        // imzML coordinates are 1-based; a <1 value is non-conformant. Warn and skip the
        // spectrum (it stays accessible by index) rather than aborting.
        OPENMS_LOG_WARN << "imaging: pixel coordinates must be >= 1; skipping spectrum " << i
                        << " at (" << StringUtils::toStr(x_imz) << "," << StringUtils::toStr(y_imz) << ")." << std::endl;
        continue;
      }
      if (spec.metaValueExists(idx.z) && static_cast<Int>(spec.getMetaValue(idx.z)) != 1)
      {
        continue; // MSImagingGeometry is 2D: only the first plane is mapped
      }

      const UInt x = static_cast<UInt>(x_imz - 1);
      const UInt y = static_cast<UInt>(y_imz - 1);
      // Soften the geometry's strict in-bounds rule to a warning: a pixel beyond the
      // declared grid is dropped from the geometry (the spectrum stays accessible by index).
      if (width > 0 && height > 0 && (x >= width || y >= height))
      {
        OPENMS_LOG_WARN << "imaging: pixel (" << StringUtils::toStr(x_imz) << "," << StringUtils::toStr(y_imz)
                        << ") at spectrum " << i << " lies outside the declared " << width << "x" << height
                        << " grid; excluding it from the imaging geometry." << std::endl;
        continue;
      }
      max_x = std::max(max_x, x);
      max_y = std::max(max_y, y);
      if (geom.hasPixel(x, y))
      {
        ++duplicate_count;
        if (duplicate_spectra.size() < max_listed)
        {
          duplicate_spectra.push_back(i);
        }
        continue;
      }
      geom.addPixel(x, y, i);
    }
    if (duplicate_count > 0)
    {
      const Size listed = duplicate_spectra.size();
      std::string indices;
      for (Size k = 0; k < listed; ++k) // concat spectrum ids and coords for warning
      {
        if (k > 0)
        {
          indices += ", ";
        }
        const MSSpectrum& d = exp[duplicate_spectra[k]];
        indices += StringUtils::toStr(duplicate_spectra[k])
                   + " (x=" + d.getMetaValue(idx.x).toString()
                   + ",y=" + d.getMetaValue(idx.y).toString() + ")";
      }
      if (listed < duplicate_count)
      {
        indices += ", ... and " + StringUtils::toStr(duplicate_count - listed) + " more";
      }
      OPENMS_LOG_WARN << "imaging: " << duplicate_count << " of " << n_spectra
                      << " spectra reuse a pixel coordinate already taken by an earlier spectrum. "
                      << "Only the first spectrum per pixel is mapped into the imaging geometry; "
                      << "all spectra remain reachable by index. Affected spectra: "
                      << indices << "." << std::endl; // std::endl: OPENMS_LOG_* only distributes a line on flush
    }

    // No declared grid: size it to the bounding box of what was mapped.
    UInt new_width = width;
    UInt new_height = height;
    if (new_width == 0 && (max_x > 0 || geom.getNumberOfPixels() > 0))
    {
      new_width = max_x + 1;
    }
    if (new_height == 0 && (max_y > 0 || geom.getNumberOfPixels() > 0))
    {
      new_height = max_y + 1;
    }
    if (new_width > 0 && new_height > 0 && (geom.getWidth() != new_width || geom.getHeight() != new_height))
    {
      geom.setDimensions(new_width, new_height);
    }
  }

  void MSImagingExperiment::rebuildGeometry()
  {
    // A pixel added through getGeometry().addPixel() left no coordinate on its spectrum, so
    // the rebuild cannot bring it back: say so instead of dropping it silently.
    const PixelMetaIndices idx = pixelMetaIndices_();
    const Size n_spectra = experiment_.getNrSpectra();
    Size unannotated = 0;
    for (const auto& p : geometry_.getPixels())
    {
      if (p.spectrum_index < n_spectra
          && (!experiment_[p.spectrum_index].metaValueExists(idx.x) || !experiment_[p.spectrum_index].metaValueExists(idx.y)))
      {
        ++unannotated;
      }
    }
    if (unannotated > 0)
    {
      OPENMS_LOG_WARN << "imaging: rebuildGeometry() drops " << unannotated << " pixel(s) bound to spectra that carry no "
                      << META_PIXEL_X << "/" << META_PIXEL_Y << " coordinate; bind pixels with bindPixel() or setGeometry() to keep them."
                      << std::endl; // std::endl: OPENMS_LOG_* only distributes a line on flush
    }
    geometry_.clearPixels();
    bindPixelsFromSpectra(experiment_, geometry_);
  }

  void MSImagingExperiment::validate() const
  {
    const PixelMetaIndices idx = pixelMetaIndices_();
    const Size n_spectra = experiment_.getNrSpectra();
    for (const auto& p : geometry_.getPixels())
    {
      if (p.spectrum_index >= n_spectra)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Pixel references missing spectrum",
                                      StringUtils::toStr(p.spectrum_index));
      }
      UInt x = 0;
      UInt y = 0;
      const MSSpectrum& spec = experiment_[p.spectrum_index];
      if (spec.metaValueExists(idx.x) && spec.metaValueExists(idx.y)
          && (!readPixelCoordinate_(spec, idx, x, y) || x != p.x || y != p.y))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Pixel (" + coordToString_(p.x, p.y) + ") is bound to spectrum "
                                        + StringUtils::toStr(p.spectrum_index)
                                        + ", which carries a different pixel coordinate; the spectrum list was "
                                          "reordered or replaced after the geometry was built (call rebuildGeometry())",
                                      spec.getMetaValue(idx.x).toString() + "," + spec.getMetaValue(idx.y).toString());
      }
    }
  }

} // namespace OpenMS
