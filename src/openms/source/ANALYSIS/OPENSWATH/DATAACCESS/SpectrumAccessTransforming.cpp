// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>

#include <utility>

namespace OpenMS
{

  SpectrumAccessTransforming::SpectrumAccessTransforming(OpenSwath::SpectrumAccessPtr sptr) :
    sptr_(std::move(sptr))
  {}

  SpectrumAccessTransforming::~SpectrumAccessTransforming() = default;

  size_t SpectrumAccessTransforming::getNrChromatograms() const
  {
    return sptr_->getNrChromatograms();
  }


  OpenSwath::SpectrumPtr SpectrumAccessTransforming::getSpectrumById(int id)
  {
    return sptr_->getSpectrumById(id);
  }

  OpenSwath::SpectrumMeta SpectrumAccessTransforming::getSpectrumMetaById(int id) const
  {
    return sptr_->getSpectrumMetaById(id);
  }

  std::vector<std::size_t> SpectrumAccessTransforming::getSpectraByRT(double RT, double deltaRT) const
  {
    return sptr_->getSpectraByRT(RT, deltaRT);
  }

  size_t SpectrumAccessTransforming::getNrSpectra() const
  {
    return sptr_->getNrSpectra();
  }

  OpenSwath::ChromatogramPtr SpectrumAccessTransforming::getChromatogramById(int id)
  {
    return sptr_->getChromatogramById(id);
  }

  std::string SpectrumAccessTransforming::getChromatogramNativeID(int id) const
  {
    return sptr_->getChromatogramNativeID(id);
  }

}

