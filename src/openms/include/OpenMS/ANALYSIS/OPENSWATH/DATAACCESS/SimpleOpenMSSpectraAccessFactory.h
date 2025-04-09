// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors:  Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#pragma once

#include <fstream>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <boost/shared_ptr.hpp>

namespace OpenMS
{
  /**
    @brief A factory method that returns two ISpectrumAccess implementations
  */
  class OPENMS_DLLAPI SimpleOpenMSSpectraFactory
  {
  public:

    /// Simple Factory method to get a SpectrumAccess Ptr from an MSExperiment
    static OpenSwath::SpectrumAccessPtr getSpectrumAccessOpenMSPtr(const boost::shared_ptr<OpenMS::PeakMap>& exp);

  private:

    static bool isExperimentCached(const boost::shared_ptr<OpenMS::PeakMap>& exp);
  };
}


