// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

  /**
    @brief NLargest removes all but the n largest peaks

    @htmlinclude OpenMS_NLargest.parameters

    @ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI NLargest :
    public DefaultParamHandler
  {
public:

    // @name Constructors and Destructors
    // @{

    /// default constructor
    NLargest();
    /// detailed constructor
    NLargest(UInt n);
    /// destructor
    ~NLargest() override;

    /// copy constructor
    NLargest(const NLargest & source);
    /// assignment operator
    NLargest & operator=(const NLargest & source);

    // @}

    ///
    template <typename SpectrumType>
    void filterSpectrum(SpectrumType & spectrum)
    {
      if (spectrum.size() <= peakcount_) return;

      // sort by reverse intensity
      spectrum.sortByIntensity(true);

      // keep the n largest peaks if more than n are present
      std::vector<Size> indices;
      for (Size i = 0; i != peakcount_; ++i)
      {
        indices.push_back(i);
      }
      spectrum.select(indices);
    }

    void filterPeakSpectrum(PeakSpectrum & spectrum);

    void filterPeakMap(PeakMap & exp);

    //TODO reimplement DefaultParamHandler::updateMembers_()

    // @}

protected:
    void updateMembers_() override;
    UInt peakcount_;

    /// handles the initialization of the default parameters for the 2 constructors
    void init_();

  };

}
