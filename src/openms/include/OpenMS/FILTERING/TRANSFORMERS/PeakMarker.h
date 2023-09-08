// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <map>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS
{
  /**
  @brief PeakMarker marks peaks that seem to fulfill some criterion

*/
  class OPENMS_DLLAPI PeakMarker :
    public DefaultParamHandler
  {
public:

    /// default constructor
    PeakMarker();

    /// copy constructor
    PeakMarker(const PeakMarker & source);

    /// destructor
    ~PeakMarker() override;

    /// assignment operator
    PeakMarker & operator=(const PeakMarker & source);

    /// method to mark peaks
    template <typename SpectrumType>
    void apply(std::map<double, bool> & /* marked */, SpectrumType & /* spectrum */) {}

    ///
    static const String getProductName()
    {
      return "PeakMarker";
    }

  };

}
