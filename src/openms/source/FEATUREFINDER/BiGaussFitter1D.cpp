// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/BiGaussFitter1D.h>

#include <OpenMS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/FEATUREFINDER/BiGaussModel.h>

namespace OpenMS
{
  BiGaussFitter1D::BiGaussFitter1D() :
    MaxLikeliFitter1D()
  {
    setName("BiGaussFitter1D");

    defaults_.setValue("statistics:variance1", 1.0, "Variance of the first gaussian, used for the lower half of the model.", {"advanced"});
    defaults_.setValue("statistics:variance2", 1.0, "Variance of the second gaussian, used for the upper half of the model.", {"advanced"});

    defaultsToParam_();
  }

  BiGaussFitter1D::BiGaussFitter1D(const BiGaussFitter1D& source) :
    MaxLikeliFitter1D(source)
  {
    updateMembers_();
  }

  BiGaussFitter1D::~BiGaussFitter1D() = default;

  BiGaussFitter1D& BiGaussFitter1D::operator=(const BiGaussFitter1D& source)
  {
    if (&source == this)
    {
      return *this;
    }
    MaxLikeliFitter1D::operator=(source);    
    updateMembers_();

    return *this;
  }

  BiGaussFitter1D::QualityType BiGaussFitter1D::fit1d(const RawDataArrayType& set, std::unique_ptr<InterpolationModel>& model)
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
    const CoordinateType stdev1 = sqrt(statistics1_.variance()) * tolerance_stdev_box_;
    const CoordinateType stdev2 = sqrt(statistics2_.variance()) * tolerance_stdev_box_;
    min_bb -= stdev1;
    max_bb += stdev2;


    // build model
    model = std::unique_ptr<BiGaussModel>(new BiGaussModel());
    model->setInterpolationStep(interpolation_step_);
    Param tmp;
    tmp.setValue("bounding_box:min", min_bb);
    tmp.setValue("bounding_box:max", max_bb);
    tmp.setValue("statistics:mean", statistics1_.mean());
    tmp.setValue("statistics:variance1", statistics1_.variance());
    tmp.setValue("statistics:variance2", statistics2_.variance());
    model->setParameters(tmp);

    // fit offset
    QualityType quality;
    quality = fitOffset_(model, set, stdev1, stdev2, interpolation_step_);
    if (std::isnan(quality))
    {
      quality = -1.0;
    }
    return quality;
  }

  void BiGaussFitter1D::updateMembers_()
  {
    MaxLikeliFitter1D::updateMembers_();
    statistics1_.setMean(param_.getValue("statistics:mean"));
    statistics1_.setVariance(param_.getValue("statistics:variance1"));
    statistics2_.setMean(param_.getValue("statistics:mean"));
    statistics2_.setVariance(param_.getValue("statistics:variance2"));
  }

}
