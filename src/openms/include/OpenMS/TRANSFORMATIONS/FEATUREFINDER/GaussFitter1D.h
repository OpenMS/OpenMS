// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MaxLikeliFitter1D.h>

namespace OpenMS
{
  /**
    @brief Gaussian distribution fitter (1-dim.) approximated using linear interpolation.

    @htmlinclude OpenMS_GaussFitter1D.parameters
  */
  class OPENMS_DLLAPI GaussFitter1D :
    public MaxLikeliFitter1D
  {
public:

    /// Default constructor
    GaussFitter1D();

    /// copy constructor
    GaussFitter1D(const GaussFitter1D & source);

    /// destructor
    ~GaussFitter1D() override;

    /// assignment operator
    virtual GaussFitter1D & operator=(const GaussFitter1D & source);

    /// create new GaussFitter1D object (function needed by Factory)
    static Fitter1D * create()
    {
      return new GaussFitter1D();
    }

    /// return interpolation model
    QualityType fit1d(const RawDataArrayType & range, std::unique_ptr<InterpolationModel>& model) override;

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
      return "GaussFitter1D";
    }

protected:

    void updateMembers_() override;
  };
}

