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
    @brief total intensity of peak pairs that could result from complementing fragments of charge state 1

        @htmlinclude OpenMS_ComplementFilter.parameters

        @ingroup SpectraFilter
  */
  class OPENMS_DLLAPI ComplementFilter :
    public FilterFunctor
  {
public:

    // @name Constructors and Destructors
    //@{
    /// standard constructor
    ComplementFilter();

    /// copy constructor
    ComplementFilter(const ComplementFilter & source);

    /// destructor
    ~ComplementFilter() override;
    //@}

    // @name Operators
    //@{
    /// assignment operator
    ComplementFilter & operator=(const ComplementFilter & source);
    //@}

    // @name Accessors
    //@{
    static FilterFunctor * create() { return new ComplementFilter(); }

    /// returns the total intensity of peak pairs which could result from complementing fragments
    template <typename SpectrumType>
    double apply(SpectrumType & spectrum)
    {
      if (spectrum.size() < 2)
      {
        return 0;
      }
      double tolerance = (double)param_.getValue("tolerance");
      double parentmass = 0.0;
      if (!spectrum.getPrecursors().empty()) parentmass = spectrum.getPrecursors()[0].getMZ();
      double result(0);

      spectrum.sortByPosition();

      /// @improvement think about an correct fast algorithm, not just an heuristic (Andreas)
      Size j = spectrum.size() - 1;
      for (Size i = 0; i < spectrum.size() && i <= j; /*++i*/)
      {
        double sum = spectrum[i].getPosition()[0] + spectrum[j].getPosition()[0];

        if (std::fabs(sum - parentmass) < tolerance)
        {
          result += spectrum[i].getIntensity() + spectrum[j].getIntensity();
        }

        if (sum < parentmass)
        {
          ++i;
        }
        else
        {
          if (sum > parentmass)
          {
            --j;
          }
        }
      }

      return result;
    }

    /// returns the name for registration at the factory
    static const String getProductName()
    {
      return "ComplementFilter";
    }

    //@}

  };
}
