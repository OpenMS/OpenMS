// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/FEATUREFINDER/BaseModel.h>
#include <OpenMS/ML/INTERPOLATION/LinearInterpolation.h>

namespace OpenMS
{
  /**
    @brief Abstract class for 1D-models that are approximated using linear interpolation

        Model wrapping LinearInterpolation for speed-up in calculation of predicted intensities
        Derived classes have to implement setSamples()

    @htmlinclude OpenMS_InterpolationModel.parameters

        @ingroup FeatureFinder

    */
  class OPENMS_DLLAPI InterpolationModel :
    public BaseModel
  {

public:
    typedef double IntensityType;
    typedef DPosition<1> PositionType;
    typedef double CoordinateType;
    using KeyType = double;
    typedef Math::LinearInterpolation<KeyType> LinearInterpolation;
    /// Default constructor
    InterpolationModel();
    /// copy constructor
    InterpolationModel(const InterpolationModel & source);
    /// destructor
    ~InterpolationModel() override;
    /// assignment operator
    InterpolationModel & operator=(const InterpolationModel & source);
    /// access model predicted intensity at position @p pos
    IntensityType getIntensity(const PositionType & pos) const override;
    /// access model predicted intensity at position @p pos
    IntensityType getIntensity(CoordinateType coord) const;
    /// Returns the interpolation class
    const LinearInterpolation & getInterpolation() const;
    CoordinateType getScalingFactor() const;
    virtual void setOffset(CoordinateType offset);
    void getSamples(SamplesType & cont) const override;
    /// "center" of the model, particular definition (depends on the derived model)
    virtual CoordinateType getCenter() const;
    /// set sample/supporting points of interpolation wrt params.
    virtual void setSamples();
    void setInterpolationStep(CoordinateType interpolation_step);
    void setScalingFactor(CoordinateType scaling);

protected:
    LinearInterpolation interpolation_;
    CoordinateType interpolation_step_;
    CoordinateType scaling_;
    void updateMembers_() override;
   

  };
}
