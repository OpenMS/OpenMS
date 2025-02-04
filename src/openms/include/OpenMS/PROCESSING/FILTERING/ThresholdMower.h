// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
  /**
    @brief Removes all peaks below an intensity threshold.

    @htmlinclude OpenMS_ThresholdMower.parameters

    @ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI ThresholdMower :
    public DefaultParamHandler
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    ThresholdMower();
    /// destructor
    ~ThresholdMower() override;

    /// copy constructor
    ThresholdMower(const ThresholdMower & source);
    /// assignment operator
    ThresholdMower & operator=(const ThresholdMower & source);
    // @}

    // @name Accessors
    // @{
    ///

    ///
    template <typename SpectrumType>
    void filterSpectrum(SpectrumType & spectrum)
    {
      threshold_ = ((double)param_.getValue("threshold"));
      std::vector<Size> indices;
      for (Size i = 0; i != spectrum.size(); ++i)
      {
        if (spectrum[i].getIntensity() >= threshold_)
        {
          indices.push_back(i);
        } 
      }
      spectrum.select(indices);
    }

    void filterPeakSpectrum(PeakSpectrum & spectrum);

    void filterPeakMap(PeakMap & exp);

    //TODO reimplement DefaultParamHandler::updateMembers_()

private:
    double threshold_;

    // @}
  };

}

