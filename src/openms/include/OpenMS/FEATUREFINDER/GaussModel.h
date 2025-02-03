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
      @brief Normal distribution approximated using linear interpolation

      @htmlinclude OpenMS_GaussModel.parameters
  */
  class OPENMS_DLLAPI GaussModel :
    public InterpolationModel
  {

public:
    typedef InterpolationModel::CoordinateType CoordinateType;
    typedef Math::BasicStatistics<CoordinateType> BasicStatistics;
    typedef InterpolationModel InterpolationModel;

    /// Default constructor
    GaussModel();

    /// copy constructor
    GaussModel(const GaussModel & source);

    /// destructor
    ~GaussModel() override;

    /// assignment operator
    virtual GaussModel & operator=(const GaussModel & source);

    /** @brief set the offset of the model

        The whole model will be shifted to the new offset without being computing all over.
        and without any discrepancy.
    */
    void setOffset(CoordinateType offset) override;

    /// set sample/supporting points of interpolation
    void setSamples() override;

    /// get the center of the Gaussian model i.e. the position of the maximum
    CoordinateType getCenter() const override;

protected:
    CoordinateType  min_;
    CoordinateType  max_;
    BasicStatistics statistics_;

    void updateMembers_() override;
  };
}

