// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/GaussFitter1D.h>

#include <OpenMS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/FEATUREFINDER/GaussModel.h>

namespace OpenMS
{
  GaussFitter1D::GaussFitter1D() :
    MaxLikeliFitter1D()
  {
    setName("GaussFitter1D");

    defaults_.setValue("statistics:variance", 1.0, "Variance of the model.", {"advanced"});
    defaults_.setValue("statistics:mean", 1.0, "Mean value of the model.", {"advanced"});
    defaultsToParam_();
  }

  GaussFitter1D::GaussFitter1D(const GaussFitter1D& source) :
    MaxLikeliFitter1D(source)
  {
    updateMembers_();
  }

  GaussFitter1D::~GaussFitter1D() = default;

  GaussFitter1D& GaussFitter1D::operator=(const GaussFitter1D& source)
  {
    if (&source == this)
    {
      return *this;
    }
    MaxLikeliFitter1D::operator=(source);
    updateMembers_();

    return *this;
  }

  GaussFitter1D::QualityType GaussFitter1D::fit1d(const RawDataArrayType& set, std::unique_ptr<InterpolationModel>& model)
  {
    // Calculate bounding box
    CoordinateType min_bb = set[0].getPos(), max_bb = set[0].getPos();
    for (UInt pos = 1; pos < set.size(); ++pos)
    {
      CoordinateType tmp = set[pos].getPos();
      if (min_bb > tmp)
      {
        min_bb = tmp;
      }
      if (max_bb < tmp)
      {
        max_bb = tmp;
      }
    }

    // Enlarge the bounding box by a few multiples of the standard deviation
    const CoordinateType stdev = sqrt(statistics_.variance()) * tolerance_stdev_box_;
    min_bb -= stdev;
    max_bb += stdev;


    // build model
    model = std::unique_ptr<InterpolationModel>(new GaussModel());
    
    model->setInterpolationStep(interpolation_step_);

    Param tmp;
    tmp.setValue("bounding_box:min", min_bb);
    tmp.setValue("bounding_box:max", max_bb);
    tmp.setValue("statistics:mean", statistics_.mean());
    tmp.setValue("statistics:variance", statistics_.variance());
    model->setParameters(tmp);

    // fit offset
    QualityType quality;
    quality = fitOffset_(model, set, stdev, stdev, interpolation_step_);
    if (std::isnan(quality))
    {
      quality = -1.0;
    }
    return quality;
  }

  void GaussFitter1D::updateMembers_()
  {
    MaxLikeliFitter1D::updateMembers_();
    statistics_.setMean(param_.getValue("statistics:mean"));
    statistics_.setVariance(param_.getValue("statistics:variance"));
  }

}
