// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <map>

namespace OpenMS
{
  /**
    @brief Scales each peak by ranking the peaks per spectrum and assigning intensity according to rank

        @ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI Scaler :
    public DefaultParamHandler
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    Scaler();
    /// destructor
    ~Scaler() override;

    /// copy constructor
    Scaler(const Scaler & source);
    /// assignment operator
    Scaler & operator=(const Scaler & source);

    // @}

    // @name Accessors
    // @{

    template <typename SpectrumType>
    void filterSpectrum(SpectrumType & spectrum)
    {
      if (spectrum.empty()) return;

      spectrum.sortByIntensity();
      typename SpectrumType::size_type count = spectrum.size();
      ++count;
      typename SpectrumType::PeakType::IntensityType last_int = 0.0;
      typename SpectrumType::Iterator it = spectrum.end();
      do
      {
        --it;
        if (it->getIntensity() != last_int)
        {
          --count;
        }
        last_int = it->getIntensity();
        it->setIntensity(count);
      }
      while (it != spectrum.begin());
    }

    void filterPeakSpectrum(PeakSpectrum & spectrum);

    void filterPeakMap(PeakMap & exp);

    //TODO reimplement DefaultParamHandler::updateMembers_() when introducing member variables

    // @}

  };

}
