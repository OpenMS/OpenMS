// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>

#include <utility>

namespace OpenMS
{
  SpectrumAccessOpenMS::SpectrumAccessOpenMS(boost::shared_ptr<MSExperimentType> ms_experiment)
  {
    // store shared pointer to the actual MSExperiment
    ms_experiment_ = std::move(ms_experiment);
  }

  SpectrumAccessOpenMS::~SpectrumAccessOpenMS() = default;

  SpectrumAccessOpenMS::SpectrumAccessOpenMS(const SpectrumAccessOpenMS & rhs) :
    ms_experiment_(rhs.ms_experiment_)
  {
    // this only copies the pointers and not the actual data ... 
  }


  boost::shared_ptr<OpenSwath::ISpectrumAccess> SpectrumAccessOpenMS::lightClone() const
  {
    return boost::shared_ptr<SpectrumAccessOpenMS>(new SpectrumAccessOpenMS(*this));
  }

  OpenSwath::SpectrumPtr SpectrumAccessOpenMS::getSpectrumById(int id)
  {
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrSpectra(), "Id cannot be larger than number of spectra");

    const MSSpectrumType& spectrum = (*ms_experiment_)[id];
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
    mz_array->data.reserve(spectrum.size());
    intensity_array->data.reserve(spectrum.size());
    for (const auto& it : spectrum)
    {
      mz_array->data.push_back(it.getMZ());
      intensity_array->data.push_back(it.getIntensity());
    }

    OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
    sptr->setMZArray(mz_array);
    sptr->setIntensityArray(intensity_array);

    for (const auto& fda : spectrum.getFloatDataArrays() )
    {
      OpenSwath::BinaryDataArrayPtr tmp(new OpenSwath::BinaryDataArray);
      tmp->data.reserve(fda.size());
      for (const auto& val : fda)
      {
        tmp->data.push_back(val);
      }
      tmp->description = fda.getName();
      sptr->getDataArrays().push_back(tmp);
    }

    for (const auto& ida : spectrum.getIntegerDataArrays() )
    {
      OpenSwath::BinaryDataArrayPtr tmp(new OpenSwath::BinaryDataArray);
      tmp->data.reserve(ida.size());
      for (const auto& val : ida)
      {
        tmp->data.push_back(val);
      }
      tmp->description = ida.getName();
      sptr->getDataArrays().push_back(tmp);
    }

    return sptr;
  }

  OpenSwath::SpectrumMeta SpectrumAccessOpenMS::getSpectrumMetaById(int id) const
  {
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrSpectra(), "Id cannot be larger than number of spectra");

    OpenSwath::SpectrumMeta meta;
    meta.RT = (*ms_experiment_)[id].getRT();
    meta.ms_level = (*ms_experiment_)[id].getMSLevel();
    return meta;
  }

  OpenSwath::ChromatogramPtr SpectrumAccessOpenMS::getChromatogramById(int id)
  {
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrChromatograms(), "Id cannot be larger than number of chromatograms");

    const MSChromatogramType& chromatogram = ms_experiment_->getChromatograms()[id];
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr rt_array(new OpenSwath::BinaryDataArray);
    rt_array->data.reserve(chromatogram.size());
    intensity_array->data.reserve(chromatogram.size());
    for (const auto& it : chromatogram)
    {
      rt_array->data.push_back(it.getRT());
      intensity_array->data.push_back(it.getIntensity());
    }

    OpenSwath::ChromatogramPtr cptr(new OpenSwath::Chromatogram);
    cptr->setTimeArray(rt_array);
    cptr->setIntensityArray(intensity_array);

    for (const auto& fda : chromatogram.getFloatDataArrays() )
    {
      OpenSwath::BinaryDataArrayPtr tmp(new OpenSwath::BinaryDataArray);
      tmp->data.reserve(fda.size());
      for (const auto& val : fda)
      {
        tmp->data.push_back(val);
      }
      tmp->description = fda.getName();
      cptr->getDataArrays().push_back(tmp);
    }

    for (const auto& ida : chromatogram.getIntegerDataArrays() )
    {
      OpenSwath::BinaryDataArrayPtr tmp(new OpenSwath::BinaryDataArray);
      tmp->data.reserve(ida.size());
      for (const auto& val : ida)
      {
        tmp->data.push_back(val);
      }
      tmp->description = ida.getName();
      cptr->getDataArrays().push_back(tmp);
    }

    return cptr;
  }

  std::vector<std::size_t> SpectrumAccessOpenMS::getSpectraByRT(double RT, double deltaRT) const
  {
    OPENMS_PRECONDITION(deltaRT >= 0, "Delta RT needs to be a positive number");

    // we first perform a search for the spectrum that is past the
    // beginning of the RT domain. Then we add this spectrum and try to add
    // further spectra as long as they are below RT + deltaRT.
    std::vector<std::size_t> result;
    auto spectrum = ms_experiment_->RTBegin(RT - deltaRT);
    if (spectrum == ms_experiment_->end()) return result;

    result.push_back(std::distance(ms_experiment_->begin(), spectrum));
    spectrum++;

    while (spectrum != ms_experiment_->end() && spectrum->getRT() <= RT + deltaRT)
    {
      result.push_back(spectrum - ms_experiment_->begin());
      spectrum++;
    }
    return result;
  }

  size_t SpectrumAccessOpenMS::getNrChromatograms() const
  {
    return ms_experiment_->getChromatograms().size();
  }

  ChromatogramSettings SpectrumAccessOpenMS::getChromatogramMetaInfo(int id) const
  {
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrChromatograms(), "Id cannot be larger than number of spectra");
    return ms_experiment_->getChromatograms()[id];
  }

  std::string SpectrumAccessOpenMS::getChromatogramNativeID(int id) const
  {
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrChromatograms(), "Id cannot be larger than number of spectra");
    return ms_experiment_->getChromatograms()[id].getNativeID();
  }

  size_t SpectrumAccessOpenMS::getNrSpectra() const
  {
    return ms_experiment_->size();
  }

  SpectrumSettings SpectrumAccessOpenMS::getSpectraMetaInfo(int id) const
  {
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrSpectra(), "Id cannot be larger than number of spectra");
    return (*ms_experiment_)[id];
  }

} //end namespace OpenMS
