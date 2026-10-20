// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OnDiscBrukerTimsExperiment.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/HANDLERS/BrukerTimsSpectrumReader.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>

#include <algorithm>
#include <unordered_map>
#include <utility>

namespace OpenMS
{
  /// Immutable after openFile(); shared by all copies
  struct OnDiscBrukerTimsExperiment::Index
  {
    std::string filename;
    Config config;
    std::shared_ptr<const Internal::BrukerTimsRun> run;
    std::shared_ptr<const ExperimentalSettings> settings;
    std::unordered_map<std::string, Size> by_native_id;
  };

  OnDiscBrukerTimsExperiment::OnDiscBrukerTimsExperiment() = default;

  OnDiscBrukerTimsExperiment::OnDiscBrukerTimsExperiment(const OnDiscBrukerTimsExperiment& source) :
    IMSExperimentRandomAccess(source),
    index_(source.index_)
  {
    if (index_) reader_ = std::make_unique<Internal::BrukerTimsSpectrumReader>(index_->run);
  }

  OnDiscBrukerTimsExperiment::OnDiscBrukerTimsExperiment(OnDiscBrukerTimsExperiment&& source) noexcept = default;

  OnDiscBrukerTimsExperiment& OnDiscBrukerTimsExperiment::operator=(const OnDiscBrukerTimsExperiment& source)
  {
    if (this != &source)
    {
      OnDiscBrukerTimsExperiment copy(source);
      *this = std::move(copy);
    }
    return *this;
  }

  OnDiscBrukerTimsExperiment& OnDiscBrukerTimsExperiment::operator=(OnDiscBrukerTimsExperiment&& source) noexcept = default;

  OnDiscBrukerTimsExperiment::~OnDiscBrukerTimsExperiment() = default;

  void OnDiscBrukerTimsExperiment::openFile(const std::string& path, const Config& config)
  {
    auto index = std::make_shared<Index>();
    index->filename = path;
    index->config = config;
    index->run = Internal::BrukerTimsRun::open(path, config);
    index->settings = std::make_shared<const ExperimentalSettings>(index->run->getSettings());
    index->by_native_id.reserve(index->run->size());
    for (Size i = 0; i < index->run->size(); ++i)
    {
      index->by_native_id.emplace(index->run->getNativeID(i), i);
    }
    reader_ = std::make_unique<Internal::BrukerTimsSpectrumReader>(index->run);
    index_ = std::move(index);
  }

  bool OnDiscBrukerTimsExperiment::isOpen() const
  {
    return index_ != nullptr;
  }

  const std::string& OnDiscBrukerTimsExperiment::getFilename() const
  {
    static const std::string empty;
    return index_ ? index_->filename : empty;
  }

  const OnDiscBrukerTimsExperiment::Config& OnDiscBrukerTimsExperiment::getConfig() const
  {
    static const Config defaults;
    return index_ ? index_->config : defaults;
  }

  Size OnDiscBrukerTimsExperiment::getNrSpectra() const
  {
    return index_ ? index_->run->size() : 0;
  }

  Size OnDiscBrukerTimsExperiment::getNrChromatograms() const
  {
    return 0;
  }

  void OnDiscBrukerTimsExperiment::checkIndex_(Size index) const
  {
    if (index >= getNrSpectra())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, static_cast<SignedSize>(index), getNrSpectra());
    }
  }

  MSSpectrum OnDiscBrukerTimsExperiment::getSpectrum(Size index)
  {
    checkIndex_(index);
    MSSpectrum spectrum;
    reader_->read(index, spectrum);
    return spectrum;
  }

  MSSpectrum OnDiscBrukerTimsExperiment::getSpectrumMetadata(Size index)
  {
    checkIndex_(index);
    MSSpectrum spectrum;
    index_->run->readMetadata(index, spectrum);
    return spectrum;
  }

  MSSpectrum OnDiscBrukerTimsExperiment::getSpectrumByNativeId(const std::string& native_id)
  {
    const Size index = getSpectrumIndexByNativeId(native_id);
    if (index == getNrSpectra())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, native_id);
    }
    return getSpectrum(index);
  }

  double OnDiscBrukerTimsExperiment::getSpectrumRT(Size index) const
  {
    checkIndex_(index);
    return index_->run->getRT(index);
  }

  UInt OnDiscBrukerTimsExperiment::getSpectrumMSLevel(Size index) const
  {
    checkIndex_(index);
    return index_->run->getMSLevel(index);
  }

  const std::string& OnDiscBrukerTimsExperiment::getSpectrumNativeID(Size index) const
  {
    checkIndex_(index);
    return index_->run->getNativeID(index);
  }

  Size OnDiscBrukerTimsExperiment::getSpectrumIndexByNativeId(const std::string& native_id) const
  {
    if (!index_) return 0;
    auto it = index_->by_native_id.find(native_id);
    return it == index_->by_native_id.end() ? getNrSpectra() : it->second;
  }

  bool OnDiscBrukerTimsExperiment::isSortedByRT() const
  {
    return true;
  }

  Size OnDiscBrukerTimsExperiment::findFirstRT(double rt) const
  {
    Size lo = 0, hi = getNrSpectra();
    while (lo < hi) // the spectra are sorted by retention time, as load() sorts them
    {
      const Size mid = lo + (hi - lo) / 2;
      if (index_->run->getRT(mid) < rt) lo = mid + 1;
      else hi = mid;
    }
    return lo;
  }

  MSChromatogram OnDiscBrukerTimsExperiment::getChromatogram(Size index)
  {
    throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, static_cast<SignedSize>(index), 0);
  }

  std::shared_ptr<const ExperimentalSettings> OnDiscBrukerTimsExperiment::getExperimentalSettings() const
  {
    if (!index_) return std::make_shared<const ExperimentalSettings>();
    return index_->settings;
  }

  std::unique_ptr<IMSExperimentRandomAccess> OnDiscBrukerTimsExperiment::clone() const
  {
    return std::make_unique<OnDiscBrukerTimsExperiment>(*this);
  }
} // namespace OpenMS
