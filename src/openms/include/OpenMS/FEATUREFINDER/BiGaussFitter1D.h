// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/FEATUREFINDER/MaxLikeliFitter1D.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>


namespace OpenMS
{
  /**
    @brief BiGaussian distribution fitter (1-dim.) approximated using linear interpolation.

    @htmlinclude OpenMS_BiGaussFitter1D.parameters
  */
  class OPENMS_DLLAPI BiGaussFitter1D :
    public MaxLikeliFitter1D
  {
public:

    /// Default constructor
    BiGaussFitter1D();

    /// copy constructor
    BiGaussFitter1D(const BiGaussFitter1D & source);

    /// destructor
    ~BiGaussFitter1D() override;

    /// assignment operator
    virtual BiGaussFitter1D & operator=(const BiGaussFitter1D & source);

    /// return interpolation model
    QualityType fit1d(const RawDataArrayType & range, std::unique_ptr<InterpolationModel>& model) override;

protected:

    /// statistics for first peak site
    Math::BasicStatistics<> statistics1_;
    /// statistics for second peak site
    Math::BasicStatistics<> statistics2_;

    void updateMembers_() override;
  };
}

