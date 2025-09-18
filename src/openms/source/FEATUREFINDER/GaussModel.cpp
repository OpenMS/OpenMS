// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/GaussModel.h>

namespace OpenMS
{
  GaussModel::GaussModel() :
    InterpolationModel(),
    statistics_()
  {
    setName("GaussModel");

    defaults_.setValue("bounding_box:min", 0.0, "Lower end of bounding box enclosing the data used to fit the model.", {"advanced"});
    defaults_.setValue("bounding_box:max", 1.0, "Upper end of bounding box enclosing the data used to fit the model.", {"advanced"});
    defaults_.setValue("statistics:mean", 0.0, "Centroid position of the model (Gaussian).", {"advanced"});
    defaults_.setValue("statistics:variance", 1.0, "The variance of the Gaussian.", {"advanced"});

    defaultsToParam_();
  }

  GaussModel::GaussModel(const GaussModel & source) :
    InterpolationModel(source)
  {
    setParameters(source.getParameters());
    updateMembers_();
  }

  GaussModel::~GaussModel() = default;

  GaussModel & GaussModel::operator=(const GaussModel & source)
  {
    if (&source == this)
    {
      return *this;
    }
    setParameters(source.getParameters());
    InterpolationModel::operator=(source);
    updateMembers_();

    return *this;
  }

  void GaussModel::setSamples()
  {
    LinearInterpolation::container_type & data = interpolation_.getData();
    data.clear();
    if (max_ == min_)
    {
      return;
    }
    data.reserve(UInt((max_ - min_) / interpolation_step_ + 1));
    CoordinateType pos = min_;
    for (UInt i = 0; pos < max_; ++i)
    {
      pos = min_ + i * interpolation_step_;
      data.push_back(statistics_.normalDensity_sqrt2pi(pos));
    }
    // scale data so that integral over distribution equals one
    // multiply sum by interpolation_step_ -> rectangular approximation of integral
    IntensityType factor = scaling_ / interpolation_step_ /
                           std::accumulate(data.begin(), data.end(), IntensityType(0));

    for (LinearInterpolation::container_type::iterator it = data.begin(); it != data.end(); ++it)
      *it *= factor;
    interpolation_.setScale(interpolation_step_);
    interpolation_.setOffset(min_);
  }

  void GaussModel::updateMembers_()
  {
    InterpolationModel::updateMembers_();

    min_ = param_.getValue("bounding_box:min");
    max_ = param_.getValue("bounding_box:max");
    statistics_.setMean(param_.getValue("statistics:mean"));
    statistics_.setVariance(param_.getValue("statistics:variance"));

    setSamples();
  }

  void GaussModel::setOffset(CoordinateType offset)
  {
    double diff = offset - getInterpolation().getOffset();
    min_ += diff;
    max_ += diff;
    statistics_.setMean(statistics_.mean() + diff);

    InterpolationModel::setOffset(offset);

    param_.setValue("bounding_box:min", min_);
    param_.setValue("bounding_box:max", max_);
    param_.setValue("statistics:mean", statistics_.mean());
  }

  GaussModel::CoordinateType GaussModel::getCenter() const
  {
    return statistics_.mean();
  }

}
