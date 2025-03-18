// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
#include <OpenMS/FEATUREFINDER/InterpolationModel.h>

using namespace OpenMS;

InterpolationModel::InterpolationModel() :
    BaseModel(),
    interpolation_()
{
    this->defaults_.setValue("interpolation_step", 0.1, "Sampling rate for the interpolation of the model function ");
    this->defaults_.setValue("intensity_scaling", 1.0, "Scaling factor used to adjust the model distribution to the intensities of the data");
    defaultsToParam_();
}

InterpolationModel::InterpolationModel(const InterpolationModel & source) :
    BaseModel(source),
    interpolation_(source.interpolation_),
    interpolation_step_(source.interpolation_step_),
    scaling_(source.scaling_)
{
    updateMembers_();
}

InterpolationModel & InterpolationModel::operator=(const InterpolationModel & source)
{
    if (&source == this) return *this;

    BaseModel::operator=(source);
    interpolation_step_ = source.interpolation_step_;
    interpolation_ = source.interpolation_;
    scaling_ = source.scaling_;

    updateMembers_();

    return *this;
}

InterpolationModel::~InterpolationModel() = default;

InterpolationModel::IntensityType InterpolationModel::getIntensity(const PositionType & pos) const
{
    return interpolation_.value(pos[0]);
}

InterpolationModel::IntensityType InterpolationModel::getIntensity(CoordinateType coord) const
{
    return interpolation_.value(coord);
}

const InterpolationModel::LinearInterpolation & InterpolationModel::getInterpolation() const
{
    return interpolation_;
}

InterpolationModel::CoordinateType InterpolationModel::getScalingFactor() const
{
    return scaling_;
}

void InterpolationModel::setOffset(CoordinateType offset)
{
    interpolation_.setOffset(offset);
}

void InterpolationModel::getSamples(SamplesType & cont) const
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

InterpolationModel::CoordinateType InterpolationModel::getCenter() const
{
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
}

void InterpolationModel::setSamples()
{
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
}

void InterpolationModel::setInterpolationStep(CoordinateType interpolation_step)
{
    interpolation_step_ = interpolation_step;
    this->param_.setValue("interpolation_step", interpolation_step_);
}

void InterpolationModel::setScalingFactor(CoordinateType scaling)
{
    scaling_ = scaling;
    this->param_.setValue("intensity_scaling", scaling_);
}

void InterpolationModel::updateMembers_()
{
    BaseModel::updateMembers_();
    interpolation_step_ = this->param_.getValue("interpolation_step");
    scaling_ = this->param_.getValue("intensity_scaling");
}
