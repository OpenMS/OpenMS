// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/QC/QCBase.h>

/**
 * @brief Number of MS spectra per MS level (SpectrumCount) as a QC metric
 */

namespace OpenMS
{
  class MSExperiment;

  class OPENMS_DLLAPI SpectrumCount : public QCBase
  {
  public:
    /// Constructor
    SpectrumCount() = default;

    /// Destructor
    virtual ~SpectrumCount() = default;

    /**
    @brief Compute number of spectra per MS level and returns them in a map

    @param exp MSExperiment containing the spectra to be counted
    @return SpectrumCount
    **/

    std::map<Size, UInt> compute(const MSExperiment& exp);

    const String& getName() const override;

    QCBase::Status requirements() const override;

  private:
    const String name_ = "SpectrumCount";
  };
} // namespace OpenMS
