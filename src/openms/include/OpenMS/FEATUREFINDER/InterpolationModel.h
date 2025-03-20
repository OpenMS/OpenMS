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
    InterpolationModel() :
      BaseModel(),
      interpolation_()
    {
      this->defaults_.setValue("interpolation_step", 0.1, "Sampling rate for the interpolation of the model function ");
      this->defaults_.setValue("intensity_scaling", 1.0, "Scaling factor used to adjust the model distribution to the intensities of the data");
      defaultsToParam_();
    }

    /// copy constructor
    InterpolationModel(const InterpolationModel & source) :
      BaseModel(source),
      interpolation_(source.interpolation_),
      interpolation_step_(source.interpolation_step_),
      scaling_(source.scaling_)
    {
      updateMembers_();
    }

    /// destructor
    ~InterpolationModel() override = default;

    /// assignment operator
    InterpolationModel & operator=(const InterpolationModel & source)
    {
      if (&source == this) return *this;

      BaseModel::operator=(source);
      interpolation_step_ = source.interpolation_step_;
      interpolation_ = source.interpolation_;
      scaling_ = source.scaling_;

      updateMembers_();

      return *this;
    }

    /// access model predicted intensity at position @p pos
    IntensityType getIntensity(const PositionType & pos) const override
    {
      return interpolation_.value(pos[0]);
    }

    /// access model predicted intensity at position @p pos
    IntensityType getIntensity(CoordinateType coord) const
    {
      return interpolation_.value(coord);
    }

    /// Returns the interpolation class
    const LinearInterpolation & getInterpolation() const
    {
      return interpolation_;
    }

    /** @brief get the scaling for the model

        A scaling factor of @p scaling means that the area under the model equals
        @p scaling. Default is 1.
    */
    CoordinateType getScalingFactor() const
    {
      return scaling_;
    }

    /** @brief set the offset of the model

        The whole model will be shifted to the new offset without being recomputed all over.
        Setting takes affect immediately.
    */
    virtual void setOffset(CoordinateType offset)
    {
      interpolation_.setOffset(offset);
    }

    /// get reasonable set of samples from the model (i.e. for printing)
    void getSamples(SamplesType & cont) const override
    {
      cont.clear();
      using PeakT = BaseModel::PeakType;
      PeakT peak;
      for (Size i = 0; i < interpolation_.getData().size(); ++i)
      {
        peak.getPosition()[0] = interpolation_.index2key((KeyType)i);
        peak.setIntensity((PeakT::IntensityType)interpolation_.getData()[i]);
        cont.push_back(peak);
      }
    }

    /// "center" of the model, particular definition (depends on the derived model)
    virtual CoordinateType getCenter() const
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /// set sample/supporting points of interpolation wrt params.
    virtual void setSamples()
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /**
        @brief Set the interpolation step for the linear interpolation of the model

        For setting to take affect, call setSamples().
    */
    void setInterpolationStep(CoordinateType interpolation_step)
    {
      interpolation_step_ = interpolation_step;
      this->param_.setValue("interpolation_step", interpolation_step_);
    }

    void setScalingFactor(CoordinateType scaling)
    {
      scaling_ = scaling;
      this->param_.setValue("intensity_scaling", scaling_);
    }

protected:
    LinearInterpolation interpolation_;
    CoordinateType interpolation_step_;
    CoordinateType scaling_;

    void updateMembers_() override
    {
      BaseModel::updateMembers_();
      interpolation_step_ = this->param_.getValue("interpolation_step");
      scaling_ = this->param_.getValue("intensity_scaling");
    }

  };
}

