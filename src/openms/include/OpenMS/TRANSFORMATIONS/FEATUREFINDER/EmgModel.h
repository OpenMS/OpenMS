// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>


namespace OpenMS
{
  /**
      @brief Exponentially modified gaussian distribution model for elution profiles.

  @htmlinclude OpenMS_EmgModel.parameters
  */
  class OPENMS_DLLAPI EmgModel :
    public InterpolationModel
  {

public:
    typedef InterpolationModel::CoordinateType CoordinateType;
    typedef Math::BasicStatistics<CoordinateType> BasicStatistics;
    typedef LinearInterpolation::container_type ContainerType;

    /// Default constructor
    EmgModel();

    /// copy constructor
    EmgModel(const EmgModel & source);

    /// destructor
    ~EmgModel() override;

    /// assignment operator
    EmgModel & operator=(const EmgModel & source);

    /// create new EmgModel object (needed by Factory)
    static BaseModel<1> * create()
    {
      return new EmgModel();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
      return "EmgModel";
    }

    /// set offset without being computing all over and without any discrepancy
    void setOffset(CoordinateType offset) override;

    /// set sample/supporting points of interpolation
    void setSamples() override;

    /// get the center of the Gaussian model i.e. the position of the maximum
    CoordinateType getCenter() const override;

protected:
    CoordinateType  min_;
    CoordinateType  max_;
    BasicStatistics statistics_;
    CoordinateType height_;
    CoordinateType width_;
    CoordinateType symmetry_;
    CoordinateType retention_;

    void updateMembers_() override;
  };
}

