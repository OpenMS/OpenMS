// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FEATUREFINDER/MaxLikeliFitter1D.h>

namespace OpenMS
{
  /**
    @brief Isotope distribution fitter (1-dim.) approximated using linear interpolation.

    @htmlinclude OpenMS_IsotopeFitter1D.parameters
   */
  class OPENMS_DLLAPI IsotopeFitter1D :
    public MaxLikeliFitter1D
  {
public:

    /// Default constructor
    IsotopeFitter1D();

    /// copy constructor
    IsotopeFitter1D(const IsotopeFitter1D & source);

    /// destructor
    ~IsotopeFitter1D() override;

    /// assignment operator
    virtual IsotopeFitter1D & operator=(const IsotopeFitter1D & source);

    /// return interpolation model
    QualityType fit1d(const RawDataArrayType & range, std::unique_ptr<InterpolationModel>& model) override;

protected:

    /// isotope charge
    CoordinateType charge_;
    /// standard derivation in isotope
    CoordinateType isotope_stdev_;
    /// maximum isotopic rank to be considered
    Int max_isotope_;

    void updateMembers_() override;
  };
}

