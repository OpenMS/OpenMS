// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MSDataStoringConsumer.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>


namespace OpenMS
{
  MSDataStoringConsumer::MSDataStoringConsumer() = default;

  void MSDataStoringConsumer::setExperimentalSettings(const ExperimentalSettings & settings)
  {
    exp_ = settings; // only override the settings, keep the data
  }

  void MSDataStoringConsumer::setExpectedSize(Size s_size, Size c_size)
  {
    exp_.reserveSpaceSpectra(s_size);
    exp_.reserveSpaceChromatograms(c_size);
  }

  void MSDataStoringConsumer::consumeSpectrum(SpectrumType & s)
  {
    exp_.addSpectrum(s);
  }

  void MSDataStoringConsumer::consumeChromatogram(ChromatogramType & c)
  {
    exp_.addChromatogram(c);
  }

  const PeakMap& MSDataStoringConsumer::getData() const
  {
    return exp_;
  }
} // namespace OpenMS
