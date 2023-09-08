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
#include <cmath>

namespace OpenMS
{
  /**
    @brief IsotopeDiffFilter returns total intensity of peak pairs that could result from isotope peaks

        @htmlinclude OpenMS_IsotopeDiffFilter.parameters

        @ingroup SpectraFilter
  */
  class OPENMS_DLLAPI IsotopeDiffFilter :
    public FilterFunctor
  {

public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    IsotopeDiffFilter();

    /// copy constructor
    IsotopeDiffFilter(const IsotopeDiffFilter & source);

    /// destructor
    ~IsotopeDiffFilter() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    IsotopeDiffFilter & operator=(const IsotopeDiffFilter & source);
    // @}

    // @name Accessors
    // @{
    ///
    static FilterFunctor * create() { return new IsotopeDiffFilter(); }

    ///
    template <typename SpectrumType>
    double apply(SpectrumType & spectrum)
    {
      double tolerance = (double)param_.getValue("tolerance");
      double isodiff = 0;

      //iterate over all peaks
      for (Size i = 0; i < spectrum.size(); ++i)
      {
        for (Size j = 1; i + j < spectrum.size(); ++j)
        {
          double pos_ij = spectrum[i + j].getPosition()[0];
          double pos_i = spectrum[i].getPosition()[0];
          if (std::fabs(pos_ij - pos_i + 1) < tolerance)
          {
            isodiff += spectrum[i].getIntensity() + spectrum[i + j].getIntensity();
          }
          else
          {
            if (std::fabs(spectrum[i + j].getPosition()[0] - spectrum[i].getPosition()[0]) > 1 + tolerance)
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
      return "IsotopeDiffFilter";
    }

    // @}

private:
  };
}
