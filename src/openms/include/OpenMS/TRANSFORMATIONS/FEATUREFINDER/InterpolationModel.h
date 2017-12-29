// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_INTERPOLATIONMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_INTERPOLATIONMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>
#include <OpenMS/MATH/MISC/LinearInterpolation.h>

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
    public BaseModel<1>
  {

public:
    typedef double IntensityType;
    typedef DPosition<1> PositionType;
    typedef double CoordinateType;
    typedef Math::LinearInterpolation<double> LinearInterpolation;

    /// Default constructor
    InterpolationModel() :
      BaseModel<1>(),
      interpolation_()
    {
      this->defaults_.setValue("interpolation_step", 0.1, "Sampling rate for the interpolation of the model function ");
      this->defaults_.setValue("intensity_scaling", 1.0, "Scaling factor used to adjust the model distribution to the intensities of the data");
    }

    /// copy constructor
    InterpolationModel(const InterpolationModel & source) :
      BaseModel<1>(source),
      interpolation_(source.interpolation_),
      interpolation_step_(source.interpolation_step_),
      scaling_(source.scaling_)
    {
    }

    /// destructor
    ~InterpolationModel() override
    {
    }

    /// assignment operator
    virtual InterpolationModel & operator=(const InterpolationModel & source)
    {
      if (&source == this) return *this;

      BaseModel<1>::operator=(source);
      interpolation_step_ = source.interpolation_step_;
      interpolation_ = source.interpolation_;
      scaling_ = source.scaling_;

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
      cont = SamplesType();
      BaseModel<1>::PeakType peak;
      for (Size i = 0; i < interpolation_.getData().size(); ++i)
      {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"
        peak.setIntensity(interpolation_.getData()[i]);
#pragma clang diagnostic pop
        peak.getPosition()[0] = interpolation_.index2key(i);
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
      BaseModel<1>::updateMembers_();
      interpolation_step_ = this->param_.getValue("interpolation_step");
      scaling_ = this->param_.getValue("intensity_scaling");
    }

  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_INTERPOLATIONMODEL_H
