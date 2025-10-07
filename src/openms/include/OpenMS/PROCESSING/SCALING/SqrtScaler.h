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
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <cmath>

namespace OpenMS
{
  /**
    @brief Scales the intensity of peaks to the sqrt

        @ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI SqrtScaler :
    public DefaultParamHandler
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    SqrtScaler();
    /// destructor
    ~SqrtScaler() override;

    /// copy constructor
    SqrtScaler(const SqrtScaler & source);
    /// assignment operator
    SqrtScaler & operator=(const SqrtScaler & source);
    // @}

    ///
    template <typename SpectrumType>
    void filterSpectrum(SpectrumType & spectrum)
    {
      bool warning = false;
      for (typename SpectrumType::Iterator it = spectrum.begin(); it != spectrum.end(); ++it)
      {
        double intens = it->getIntensity();
        if (intens < 0)
        {
          intens = 0;
          warning = true;
        }
        it->setIntensity(std::sqrt(intens));
      }
      if (warning)
      {
        std::cerr << "Warning negative intensities were set to zero" << std::endl;
      }
      return;
    }

    void filterPeakSpectrum(PeakSpectrum & spectrum);

    void filterPeakMap(PeakMap & exp);

    //TODO reimplement DefaultParamHandler::updateMembers_()

    // @}

  };

}

