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
#include <OpenMS/CONCEPT/Exception.h>

#include <cstdio>
#include <functional>
#include <unordered_map>
#include <utility>

namespace OpenMS
{

struct OnDiscImzMLExperiment::Impl
{
  FILE* ibd_ {nullptr};
  ImzMLMeta meta_;
  std::vector<ImzMLSpectrumIndex> index_;

  struct PixelKey
  {
    uint32_t x {0};
    uint32_t y {0};
    uint32_t z {0};

    bool operator==(const PixelKey& other) const noexcept
    {
      return x == other.x && y == other.y && z == other.z;
    }
  };

  struct PixelKeyHash
  {
    std::size_t operator()(const PixelKey& key) const noexcept
    {
      const std::size_t hx = std::hash<uint32_t>()(key.x);
      const std::size_t hy = std::hash<uint32_t>()(key.y);
      const std::size_t hz = std::hash<uint32_t>()(key.z);
      return hx ^ (hy << 1) ^ (hz << 2);
    }
  };

  mutable std::unordered_map<PixelKey, std::size_t, PixelKeyHash> coord_map_;
  mutable bool coord_map_built_ {false};

  bool is_open_ {false};
  String ibd_path_;

  ~Impl()
  {
    if (ibd_)
    {
      fclose(ibd_);
      ibd_ = nullptr;
    }
  }

  void buildCoordMap() const
  {
    coord_map_.clear();
    coord_map_.reserve(index_.size());
    for (std::size_t i = 0; i < index_.size(); ++i)
    {
      const auto& e = index_[i];
      const PixelKey key {e.x, e.y, e.z};
      if (!coord_map_.emplace(key, i).second)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Duplicate pixel coordinate",
                                      String("(") + e.x + "," + e.y + "," + e.z + ")");
      }
    }
    coord_map_built_ = true;
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
                                  String("m/z and intensity array length mismatch at pixel (")
                                    + e.x + "," + e.y + "," + e.z + "): mz=" + mz_vec.size()
                                    + " intensity=" + int_vec.size());
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

void OnDiscImzMLExperiment::open(const String& imzml_path, const String& ibd_path_override)
{
  pimpl_ = std::make_unique<Impl>();

  ImzMLFile f;
  f.setLogType(ProgressLogger::NONE);
  f.loadSpectraIndex(imzml_path, pimpl_->meta_, pimpl_->index_);

  pimpl_->ibd_path_ = ibd_path_override.empty() ? pimpl_->meta_.ibd_file_path : ibd_path_override;
  if (pimpl_->ibd_path_.empty())
  {
    String lower = imzml_path;
    lower.toLower();
    pimpl_->ibd_path_ = lower.hasSuffix(".imzml")
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
  pimpl_->coord_map_built_ = false;
  pimpl_->coord_map_.clear();
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
  if (!pimpl_->coord_map_built_) pimpl_->buildCoordMap();

  const Impl::PixelKey key {x, y, z};
  auto it = pimpl_->coord_map_.find(key);
  if (it == pimpl_->coord_map_.end())
    throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     String("(") + x + "," + y + "," + z + ")");
  return getSpectrum(it->second);
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
