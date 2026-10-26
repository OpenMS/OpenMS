// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOnDisc.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <algorithm>
#include <numeric>

namespace OpenMS
{
  namespace
  {
    template <typename Array>
    void appendArrays(const std::vector<Array>& arrays, std::vector<OpenSwath::BinaryDataArrayPtr>& out)
    {
      for (const auto& array : arrays)
      {
        auto converted = std::make_shared<OpenSwath::BinaryDataArray>();
        converted->data.assign(array.begin(), array.end());
        converted->description = array.getName();
        out.push_back(std::move(converted));
      }
    }
  } // namespace

  SpectrumAccessOnDisc::SpectrumAccessOnDisc(std::shared_ptr<IMSExperimentRandomAccess> experiment) :
    SpectrumAccessOnDisc(experiment, [&] {
      std::vector<Size> all(experiment ? experiment->getNrSpectra() : 0);
      std::iota(all.begin(), all.end(), Size(0));
      return all;
    }())
  {
  }

  SpectrumAccessOnDisc::SpectrumAccessOnDisc(std::shared_ptr<IMSExperimentRandomAccess> experiment, std::vector<Size> indices) :
    experiment_(std::move(experiment))
  {
    if (! experiment_) { throw Exception::NullPointer(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
    auto view = std::make_shared<View>();
    const Size n = experiment_->getNrSpectra();
    std::vector<std::pair<double, Size>> by_rt;
    by_rt.reserve(indices.size());
    for (Size index : indices)
    {
      if (index >= n) { throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, static_cast<SignedSize>(index), n); }
      by_rt.emplace_back(experiment_->getSpectrumRT(index), index);
    }
    std::stable_sort(by_rt.begin(), by_rt.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    view->indices.reserve(by_rt.size());
    view->rts.reserve(by_rt.size());
    for (const auto& [rt, index] : by_rt)
    {
      view->rts.push_back(rt);
      view->indices.push_back(index);
    }
    view_ = std::move(view);
  }

  SpectrumAccessOnDisc::SpectrumAccessOnDisc(const SpectrumAccessOnDisc& rhs) :
    OpenSwath::ISpectrumAccess(rhs),
    experiment_(rhs.experiment_->clone()),
    view_(rhs.view_)
  {
  }

  SpectrumAccessOnDisc::~SpectrumAccessOnDisc() = default;

  std::shared_ptr<OpenSwath::ISpectrumAccess> SpectrumAccessOnDisc::lightClone() const
  {
    return std::make_shared<SpectrumAccessOnDisc>(*this);
  }

  Size SpectrumAccessOnDisc::getExperimentIndex(int id) const
  {
    if (id < 0 || static_cast<Size>(id) >= view_->indices.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, id, view_->indices.size());
    }
    return view_->indices[id];
  }

  const IMSExperimentRandomAccess& SpectrumAccessOnDisc::getExperiment() const
  {
    return *experiment_;
  }

  OpenSwath::SpectrumPtr SpectrumAccessOnDisc::getSpectrumById(int id)
  {
    const MSSpectrum spectrum = experiment_->getSpectrumPeaks(getExperimentIndex(id));
    auto sptr = std::make_shared<OpenSwath::Spectrum>();
    OpenSwath::BinaryDataArrayPtr mz_array = sptr->getMZArray();
    OpenSwath::BinaryDataArrayPtr intensity_array = sptr->getIntensityArray();
    mz_array->data.reserve(spectrum.size());
    intensity_array->data.reserve(spectrum.size());
    for (const auto& peak : spectrum)
    {
      mz_array->data.push_back(peak.getMZ());
      intensity_array->data.push_back(peak.getIntensity());
    }
    appendArrays(spectrum.getFloatDataArrays(), sptr->getDataArrays());
    appendArrays(spectrum.getIntegerDataArrays(), sptr->getDataArrays());
    return sptr;
  }

  OpenSwath::SpectrumMeta SpectrumAccessOnDisc::getSpectrumMetaById(int id) const
  {
    const Size index = getExperimentIndex(id);
    OpenSwath::SpectrumMeta meta;
    meta.index = index;
    meta.RT = view_->rts[id];
    meta.ms_level = static_cast<int>(experiment_->getSpectrumMSLevel(index));
    return meta;
  }

  std::vector<std::size_t> SpectrumAccessOnDisc::getSpectraByRT(double RT, double deltaRT) const
  {
    OPENMS_PRECONDITION(deltaRT >= 0, "Delta RT needs to be a positive number");
    const auto& rts = view_->rts;
    std::vector<std::size_t> result;
    auto it = std::lower_bound(rts.begin(), rts.end(), RT - deltaRT);
    if (it == rts.end()) { return result; }
    result.push_back(std::distance(rts.begin(), it));
    for (++it; it != rts.end() && *it <= RT + deltaRT; ++it)
    {
      result.push_back(std::distance(rts.begin(), it));
    }
    return result;
  }

  size_t SpectrumAccessOnDisc::getNrSpectra() const
  {
    return view_->indices.size();
  }

  OpenSwath::ChromatogramPtr SpectrumAccessOnDisc::getChromatogramById(int id)
  {
    const MSChromatogram chromatogram = experiment_->getChromatogram(static_cast<Size>(id));
    auto cptr = std::make_shared<OpenSwath::Chromatogram>();
    OpenSwath::BinaryDataArrayPtr rt_array = cptr->getTimeArray();
    OpenSwath::BinaryDataArrayPtr intensity_array = cptr->getIntensityArray();
    rt_array->data.reserve(chromatogram.size());
    intensity_array->data.reserve(chromatogram.size());
    for (const auto& peak : chromatogram)
    {
      rt_array->data.push_back(peak.getRT());
      intensity_array->data.push_back(peak.getIntensity());
    }
    appendArrays(chromatogram.getFloatDataArrays(), cptr->getDataArrays());
    appendArrays(chromatogram.getIntegerDataArrays(), cptr->getDataArrays());
    return cptr;
  }

  size_t SpectrumAccessOnDisc::getNrChromatograms() const
  {
    return experiment_->getNrChromatograms();
  }

  std::string SpectrumAccessOnDisc::getChromatogramNativeID(int id) const
  {
    return experiment_->getChromatogram(static_cast<Size>(id)).getNativeID();
  }
} // namespace OpenMS
