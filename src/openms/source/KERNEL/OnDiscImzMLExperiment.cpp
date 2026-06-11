// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Aditya Sarna $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/OnDiscImzMLExperiment.h>
#include <OpenMS/FORMAT/ImzMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/ImzMLHandlerHelper.h>
#include <OpenMS/IMAGING/MSImagingGeometry.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <algorithm>
#include <cstdio>
#include <utility>

namespace OpenMS
{

struct OnDiscImzMLExperiment::Impl
{
  FILE* ibd_ {nullptr};
  ImzMLMeta meta_;
  std::vector<ImzMLSpectrumIndex> index_;

  // Standard 2D imaging geometry (pixel grid + (x,y)->spectrum index), shared with
  // the in-memory path. Built eagerly in open(): it is an in-memory O(n) pass over
  // the already-parsed index (no .ibd reads), negligible next to the XML parse open()
  // already performs — only peak decoding is worth deferring. Building it at open()
  // also fails fast on a structurally broken coordinate grid (duplicate / <1 coords).
  MSImagingGeometry geometry_;

  bool is_open_ {false};
  std::string ibd_path_;

  ~Impl()
  {
    if (ibd_)
    {
      fclose(ibd_);
      ibd_ = nullptr;
    }
  }

  // Populate geometry_ from the parsed index. Mirrors ImzMLFile::buildImagingGeometry
  // (the in-memory sibling) so the two paths normalize coordinates identically:
  // imzML 1-based -> 0-based, z == 1 plane only, dimensions and pixel size from
  // ImzMLMeta. Keep the two implementations in sync. <1 and out-of-grid coordinates are
  // warned and skipped; only duplicate coordinates raise Exception::InvalidValue.
  void buildGeometry_()
  {
    geometry_.clear();

    UInt width = meta_.max_count_x;
    UInt height = meta_.max_count_y;
    if (width > 0 && height > 0)
    {
      geometry_.setDimensions(width, height);
    }

    UInt max_x = 0;
    UInt max_y = 0;
    for (std::size_t i = 0; i < index_.size(); ++i)
    {
      const ImzMLSpectrumIndex& e = index_[i];
      if (e.z != 1) // MSImagingGeometry is 2D: only the first plane is mapped
      {
        continue;
      }
      if (e.x < 1 || e.y < 1)
      {
        // imzML coordinates are 1-based; a <1 value is non-conformant. Warn and skip the
        // spectrum (it stays reachable by linear index) rather than aborting the load.
        OPENMS_LOG_WARN << "imzML: pixel coordinates must be >= 1; skipping spectrum " << i
                        << " at (" << StringConversions::toString(e.x) << ","
                        << StringConversions::toString(e.y) << ")." << std::endl;
        continue;
      }
      const UInt x = static_cast<UInt>(e.x - 1);
      const UInt y = static_cast<UInt>(e.y - 1);
      // Soften the geometry's strict in-bounds rule to a warning: a pixel beyond the
      // declared grid (IMS:1000042/043) is dropped from the geometry (the spectrum stays
      // reachable by linear index) instead of aborting the load.
      if (width > 0 && height > 0 && (x >= width || y >= height))
      {
        OPENMS_LOG_WARN << "imzML: pixel (" << StringConversions::toString(e.x) << ","
                        << StringConversions::toString(e.y) << ") at spectrum " << i
                        << " lies outside the declared " << width << "x" << height
                        << " grid; excluding it from the imaging geometry." << std::endl;
        continue;
      }
      max_x = std::max(max_x, x);
      max_y = std::max(max_y, y);
      geometry_.addPixel(x, y, i); // still throws Exception::InvalidValue on duplicate coordinates
    }

    if (width == 0 && (max_x > 0 || geometry_.getNumberOfPixels() > 0))
    {
      width = max_x + 1;
    }
    if (height == 0 && (max_y > 0 || geometry_.getNumberOfPixels() > 0))
    {
      height = max_y + 1;
    }
    if (width > 0 && height > 0 && (geometry_.getWidth() != width || geometry_.getHeight() != height))
    {
      geometry_.setDimensions(width, height);
    }

    if (meta_.pixel_size_x > 0 && meta_.pixel_size_y > 0)
    {
      geometry_.setPixelSize(meta_.pixel_size_x, meta_.pixel_size_y, "micrometer");
    }
  }

  MSSpectrum decodeSpectrum(std::size_t i) const
  {
    const ImzMLSpectrumIndex& e = index_[i];
    MSSpectrum s;
    s.setMetaValue("imzml:x", static_cast<Int>(e.x));
    s.setMetaValue("imzml:y", static_cast<Int>(e.y));
    s.setMetaValue("imzml:z", static_cast<Int>(e.z));

    std::vector<double> mz_vec = readMz_(e);
    std::vector<float> int_vec = readInt_(e);

    if (mz_vec.size() != int_vec.size())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, ibd_path_,
                                  "m/z and intensity array length mismatch at pixel ("
                                    + StringConversions::toString(e.x) + "," + StringConversions::toString(e.y) + "," + StringConversions::toString(e.z) + "): mz=" + StringConversions::toString(mz_vec.size())
                                    + " intensity=" + StringConversions::toString(int_vec.size()));
    }

    const std::size_t n = mz_vec.size();
    if (n == 0 && (e.mz_length > 0 || e.int_length > 0))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, ibd_path_,
                                  "failed to decode spectrum peaks from .ibd");
    }
    s.resize(n);
    for (std::size_t k = 0; k < n; ++k)
    {
      s[k].setMZ(mz_vec[k]);
      s[k].setIntensity(int_vec[k]);
    }
    s.sortByPosition();
    return s;
  }

private:
  std::vector<double> readMz_(const ImzMLSpectrumIndex& e) const
  {
    std::vector<double> out;
    if (!e.mz_length)
    {
      return out;
    }
    ImzMLBinaryIO::readMzArray(ibd_, e.mz_offset, e.mz_length, e.mz_type, out, ibd_path_);
    return out;
  }

  std::vector<float> readInt_(const ImzMLSpectrumIndex& e) const
  {
    std::vector<float> out;
    if (!e.int_length)
    {
      return out;
    }
    ImzMLBinaryIO::readIntArray(ibd_, e.int_offset, e.int_length, e.int_type, out, ibd_path_);
    return out;
  }
};

OnDiscImzMLExperiment::OnDiscImzMLExperiment()
  : pimpl_(std::make_unique<Impl>())
{
}

OnDiscImzMLExperiment::~OnDiscImzMLExperiment() = default;

OnDiscImzMLExperiment::OnDiscImzMLExperiment(OnDiscImzMLExperiment&& other)
  : pimpl_(std::move(other.pimpl_))
{
  other.pimpl_ = std::make_unique<Impl>();
}

OnDiscImzMLExperiment& OnDiscImzMLExperiment::operator=(OnDiscImzMLExperiment&& other)
{
  if (this != &other)
  {
    close();
    pimpl_ = std::move(other.pimpl_);
    other.pimpl_ = std::make_unique<Impl>();
  }
  return *this;
}

void OnDiscImzMLExperiment::open(const std::string& imzml_path, const std::string& ibd_path_override)
{
  pimpl_ = std::make_unique<Impl>();

  ImzMLFile f;
  f.setLogType(ProgressLogger::NONE);
  f.loadSpectraIndex(imzml_path, pimpl_->meta_, pimpl_->index_, ibd_path_override);

  // Build the 2D imaging geometry up front (in-memory pass over the index, no .ibd
  // reads). Fails fast here on a broken coordinate grid; only peak decode stays lazy.
  pimpl_->buildGeometry_();

  pimpl_->ibd_path_ = ibd_path_override.empty() ? pimpl_->meta_.ibd_file_path : ibd_path_override;
  if (pimpl_->ibd_path_.empty())
  {
    std::string lower = imzml_path;
    StringUtils::toLower(lower);
    pimpl_->ibd_path_ = StringUtils::hasSuffix(lower, ".imzml")
                          ? imzml_path.substr(0, imzml_path.size() - 6) + ".ibd"
                          : imzml_path + ".ibd";
  }

  pimpl_->ibd_ = fopen(pimpl_->ibd_path_.c_str(), "rb");
  if (!pimpl_->ibd_)
    throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pimpl_->ibd_path_);

  pimpl_->is_open_ = true;
}

void OnDiscImzMLExperiment::close() noexcept
{
  if (pimpl_->ibd_)
  {
    fclose(pimpl_->ibd_);
    pimpl_->ibd_ = nullptr;
  }
  pimpl_->is_open_ = false;
  pimpl_->geometry_.clear();
}

bool OnDiscImzMLExperiment::isOpen() const noexcept
{
  return pimpl_->is_open_;
}

std::size_t OnDiscImzMLExperiment::getNrSpectra() const noexcept
{
  return pimpl_->index_.size();
}

const ImzMLSpectrumIndex& OnDiscImzMLExperiment::getIndex(std::size_t i) const
{
  if (i >= pimpl_->index_.size())
    throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                   static_cast<SignedSize>(i),
                                   static_cast<SignedSize>(pimpl_->index_.size()));
  return pimpl_->index_[i];
}

MSSpectrum OnDiscImzMLExperiment::getSpectrum(std::size_t i) const
{
  if (i >= pimpl_->index_.size())
    throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                   static_cast<SignedSize>(i),
                                   static_cast<SignedSize>(pimpl_->index_.size()));
  if (!pimpl_->ibd_)
    throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pimpl_->ibd_path_);
  return pimpl_->decodeSpectrum(i);
}

MSSpectrum OnDiscImzMLExperiment::getSpectrumAtCoord(uint32_t x, uint32_t y, uint32_t z) const
{
  if (!pimpl_->ibd_)
    throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pimpl_->ibd_path_);

  // imzML coordinates are 1-based; the shared MSImagingGeometry is 0-based and models
  // only the z == 1 plane. Convert and look up there.
  if (z == 1 && x >= 1 && y >= 1)
  {
    const MSImagingGeometry& geom = pimpl_->geometry_; // built during open()
    if (geom.hasPixel(x - 1, y - 1))
    {
      return getSpectrum(geom.getSpectrumIndex(x - 1, y - 1));
    }
  }
  throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                   "(" + StringConversions::toString(x) + "," + StringConversions::toString(y) + "," + StringConversions::toString(z) + ")");
}

const MSImagingGeometry& OnDiscImzMLExperiment::getGeometry() const
{
  return pimpl_->geometry_;
}

const ImzMLMeta& OnDiscImzMLExperiment::getImzMLMeta() const noexcept
{
  return pimpl_->meta_;
}

uint32_t OnDiscImzMLExperiment::gridWidth() const noexcept
{
  return pimpl_->meta_.max_count_x;
}

uint32_t OnDiscImzMLExperiment::gridHeight() const noexcept
{
  return pimpl_->meta_.max_count_y;
}

} // namespace OpenMS
