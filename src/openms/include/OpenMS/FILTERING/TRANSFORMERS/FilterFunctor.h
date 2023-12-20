// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS
{
  /**
    @brief A FilterFunctor extracts some spectrum characteristics for quality assessment
  */
  class OPENMS_DLLAPI FilterFunctor :
    public DefaultParamHandler
  {
public:

    /// default constructor
    FilterFunctor();

    /// copy constructor
    FilterFunctor(const FilterFunctor & source);

    /// destructor
    ~FilterFunctor() override;

    /// assignment operator
    FilterFunctor & operator=(const FilterFunctor & source);

    ///
    static void registerChildren();

    /// function call operator
    template <typename SpectrumType>
    double apply(SpectrumType & /* spectrum */)
    {
      return 0;
    }

  };
}
