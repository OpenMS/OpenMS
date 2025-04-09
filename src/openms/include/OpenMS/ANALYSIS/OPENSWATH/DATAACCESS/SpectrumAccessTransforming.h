// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

namespace OpenMS
{

  /**
   * @brief An abstract base class implementing a transforming wrapper around spectrum access.
   *
   */
  class OPENMS_DLLAPI SpectrumAccessTransforming :
    public OpenSwath::ISpectrumAccess
  {
public:

    explicit SpectrumAccessTransforming(OpenSwath::SpectrumAccessPtr sptr);
        
    ~SpectrumAccessTransforming() override = 0;

    boost::shared_ptr<ISpectrumAccess> lightClone() const override = 0;

    OpenSwath::SpectrumPtr getSpectrumById(int id) override;

    OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const override;

    std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const override;

    size_t getNrSpectra() const override;

    OpenSwath::ChromatogramPtr getChromatogramById(int id) override;

    size_t getNrChromatograms() const override;

    std::string getChromatogramNativeID(int id) const override;

protected:
    OpenSwath::SpectrumAccessPtr sptr_;

  };

}


