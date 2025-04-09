// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

namespace OpenMS
{

  /**
    @brief Consumer class that performs no operation.

    This is sometimes necessary to fulfill the requirement of passing an
    valid Interfaces::IMSDataConsumer object or pointer but no operation is
    required.

  */
  class OPENMS_DLLAPI NoopMSDataConsumer :
    public Interfaces::IMSDataConsumer
  {
  public:

    NoopMSDataConsumer() {}
    void setExperimentalSettings(const ExperimentalSettings &) override {}
    void setExpectedSize(Size, Size) override {}
    void consumeSpectrum(SpectrumType &) override {}
    void consumeChromatogram(ChromatogramType &) override {}
  };
} //end namespace OpenMS


