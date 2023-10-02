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

namespace OpenMS
{
  /**
    @brief TICFilter calculates TIC

        @htmlinclude OpenMS_TICFilter.parameters

        @ingroup SpectraFilter
  */
  class OPENMS_DLLAPI TICFilter :
    public FilterFunctor
  {
public:

    // @name Constructors and Destructors
    // @{
    /// standard constructor
    TICFilter();

    /// copy constructor
    TICFilter(const TICFilter & source);

    /// destructor
    ~TICFilter() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    TICFilter & operator=(const TICFilter & source);
    // @}

    // @name Accessors
    // @{
    ///
    static FilterFunctor * create() { return new TICFilter(); }

    ///
    template <typename SpectrumType>
    double apply(SpectrumType & spectrum)
    {
      typedef typename SpectrumType::ConstIterator ConstIterator;
      double TIC = 0;
      //double window = (double)param_.getValue("window");

      for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
      {
        TIC += it->getIntensity();
      }
      return TIC;
    }

    ///
    static const String getProductName()
    {
      return "TICFilter";
    }

    // @}

  };
}
