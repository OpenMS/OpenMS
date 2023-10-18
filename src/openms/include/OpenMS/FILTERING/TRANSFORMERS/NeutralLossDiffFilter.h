// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>
#include <cmath>

namespace OpenMS
{
  /**
    @brief NeutralLossDiffFilter returns the total intensity ob peak pairs whose m/z difference can be explained by a neutral loss

        @htmlinclude OpenMS_NeutralLossDiffFilter.parameters

        @ingroup SpectraFilter
  */
  class OPENMS_DLLAPI NeutralLossDiffFilter :
    public FilterFunctor
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    NeutralLossDiffFilter();

    /// copy constructor
    NeutralLossDiffFilter(const NeutralLossDiffFilter & source);

    /// destructor
    ~NeutralLossDiffFilter() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    NeutralLossDiffFilter & operator=(const NeutralLossDiffFilter & source);
    // @}

    // @name Accessors
    // @{
    ///
    static FilterFunctor * create() { return new NeutralLossDiffFilter(); }

    ///
    template <typename SpectrumType>
    double apply(SpectrumType & spectrum)
    {
      double tolerance = (double)param_.getValue("tolerance");
      double isodiff = 0;
      //iterate over all peaks
      for (int i = 0; i < (int)spectrum.size(); ++i)
      {
        for (int j = 1; i - j >= 0; ++j)
        {
          double pos_diff = std::fabs(spectrum[i - j].getPosition()[0] - spectrum[i].getPosition()[0]);
          if (std::fabs(pos_diff - 18) < tolerance || std::fabs(pos_diff - 17) < tolerance)   // water and ammonium
          {
            isodiff += spectrum[i - j].getIntensity() + spectrum[i].getIntensity();
          }
          else
          {
            if (pos_diff > 18 + tolerance)
            {
              break;
            }
          }
        }
      }

      return isodiff;
    }

    ///
    static const String getProductName()
    {
      return "NeutralLossDiffFilter";
    }

    // @}

  };
}
