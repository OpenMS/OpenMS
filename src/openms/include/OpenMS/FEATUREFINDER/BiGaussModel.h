// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

namespace OpenMS
{
  /**
        @brief BiGaussian distribution approximated using linear interpolation.

        Asymmetric distribution realized via two normal distributions with
        different variances combined at the mean.

    @htmlinclude OpenMS_BiGaussModel.parameters
    */
  class OPENMS_DLLAPI BiGaussModel :
    public InterpolationModel
  {
public:
    typedef InterpolationModel::CoordinateType CoordinateType;

    /// Default constructor
    BiGaussModel();

    /// copy constructor
    BiGaussModel(const BiGaussModel & source);

    /// destructor
    ~BiGaussModel() override;

    /// assignment operator
    virtual BiGaussModel & operator=(const BiGaussModel & source);

    /** @brief set the offset of the model

        The whole model will be shifted to the new offset without being computing all over
        and without any discrepancy.
    */
    void setOffset(CoordinateType offset) override;

    /// set sample/supporting points of interpolation
    void setSamples() override;

    /// get the center of the BiGaussian model i.e. the position of the maximum
    CoordinateType getCenter() const override;

protected:
    CoordinateType min_;
    CoordinateType max_;
    Math::BasicStatistics<> statistics1_;
    Math::BasicStatistics<> statistics2_;

    void updateMembers_() override;
  };
}

