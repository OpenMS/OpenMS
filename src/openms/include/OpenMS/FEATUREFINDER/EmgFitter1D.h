// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FEATUREFINDER/LevMarqFitter1D.h>
#include <OpenMS/CONCEPT/Constants.h>

namespace OpenMS
{
  /**
    @brief Exponentially modified gaussian distribution fitter (1-dim.) using Levenberg-Marquardt algorithm (Eigen implementation) for parameter optimization.

    @htmlinclude OpenMS_EmgFitter1D.parameters
  */
  class OPENMS_DLLAPI EmgFitter1D :
    public LevMarqFitter1D
  {
public:

    /// Default constructor
    EmgFitter1D();

    /// copy constructor
    EmgFitter1D(const EmgFitter1D& source);

    /// destructor
    ~EmgFitter1D() override;

    /// assignment operator
    virtual EmgFitter1D& operator=(const EmgFitter1D& source);

    /// return interpolation model
    QualityType fit1d(const RawDataArrayType& range, std::unique_ptr<InterpolationModel>& model) override;

protected:
    void setInitialParameters_(const RawDataArrayType& set);
    void setInitialParametersMOM_(const RawDataArrayType& set);
    void updateMembers_() override;

    // Inner classes and member variables
    struct Data;
    class EgmFitterFunctor;
    CoordinateType height_;
    CoordinateType width_;
    CoordinateType symmetry_;
    CoordinateType retention_;
  };
}

