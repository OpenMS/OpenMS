// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

  /**
    @brief Consumer class that simply stores the data.

    This class is able to keep spectra and chromatograms passed to it in memory
    and the data can be accessed through getData()

  */
  class OPENMS_DLLAPI MSDataStoringConsumer :
    public Interfaces::IMSDataConsumer
  {
  private:
    PeakMap exp_;

  public:

    MSDataStoringConsumer();

    void setExperimentalSettings(const ExperimentalSettings & settings) override;

    void setExpectedSize(Size s_size, Size c_size) override;

    void consumeSpectrum(SpectrumType& s) override;

    void consumeChromatogram(ChromatogramType& c) override;

    const PeakMap& getData() const;

  };
} //end namespace OpenMS


