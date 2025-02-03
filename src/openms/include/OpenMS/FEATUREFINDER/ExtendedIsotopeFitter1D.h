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
  @brief Extended isotope distribution fitter (1-dim.) approximated using linear interpolation.

  @htmlinclude OpenMS_ExtendedIsotopeFitter1D.parameters
  */
  class OPENMS_DLLAPI ExtendedIsotopeFitter1D :
    public MaxLikeliFitter1D
  {
public:

    /// Default constructor
    ExtendedIsotopeFitter1D();

    /// copy constructor
    ExtendedIsotopeFitter1D(const ExtendedIsotopeFitter1D & source);

    /// destructor
    ~ExtendedIsotopeFitter1D() override;

    /// assignment operator
    virtual ExtendedIsotopeFitter1D & operator=(const ExtendedIsotopeFitter1D & source);

    /// return interpolation model
    QualityType fit1d(const RawDataArrayType & range, std::unique_ptr<InterpolationModel>& model) override;

protected:

    /// isotope charge
    CoordinateType charge_;
    /// standard derivation in isotope
    CoordinateType isotope_stdev_;
    /// monoisotopic mass
    CoordinateType monoisotopic_mz_;
    /// maximum isotopic rank to be considered
    Int max_isotope_;

    void updateMembers_() override;
  };
}

