// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <boost/shared_ptr.hpp>

namespace OpenMS
{
  /**
    @brief An implementation of the OpenSWATH Spectrum Access interface using OpenMS

  */
  class OPENMS_DLLAPI SpectrumAccessOpenMS :
    public OpenSwath::ISpectrumAccess
  {
public:
    typedef OpenMS::PeakMap MSExperimentType;
    typedef OpenMS::MSSpectrum MSSpectrumType;
    typedef OpenMS::MSChromatogram MSChromatogramType;

    /// Constructor
    explicit SpectrumAccessOpenMS(boost::shared_ptr<MSExperimentType> ms_experiment);

    /// Destructor
    ~SpectrumAccessOpenMS() override;

    /**
      @brief Copy constructor

      Performs a light copy operation when another SpectrumAccessOpenMS
      instance is given: only a copy of the pointer to the underlying
      MSExperiment is stored, so after this, both instances (rhs and *this)
      will point to the same MSExperiment.

    */
    SpectrumAccessOpenMS(const SpectrumAccessOpenMS & rhs);

    /**
      @brief Light clone operator (actual data will not get copied)

      Creates a light clone of the current instance, with the clone pointing to
      the same underlying MSExperiment.

    */
    boost::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const override;

    OpenSwath::SpectrumPtr getSpectrumById(int id) override;

    OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const override;

    std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const override;

    size_t getNrSpectra() const override;

    SpectrumSettings getSpectraMetaInfo(int id) const;

    OpenSwath::ChromatogramPtr getChromatogramById(int id) override;

    // FEATURE ?
    // ChromatogramPtr getChromatogramByPrecursorMZ(double mz, double deltaMZ);

    size_t getNrChromatograms() const override;

    ChromatogramSettings getChromatogramMetaInfo(int id) const;

    std::string getChromatogramNativeID(int id) const override;

private:
    boost::shared_ptr<MSExperimentType> ms_experiment_;

  };
} //end namespace OpenMS

