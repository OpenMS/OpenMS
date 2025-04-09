// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/IsotopeFitter1D.h>
#include <OpenMS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/FEATUREFINDER/GaussModel.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

namespace OpenMS
{

  IsotopeFitter1D::IsotopeFitter1D() :
    MaxLikeliFitter1D()
  {
    setName("IsotopeFitter1D");

    defaults_.setValue("statistics:variance", 1.0, "Variance of the model.", {"advanced"});
    defaults_.setValue("charge", 1, "Charge state of the model.", {"advanced"});
    defaults_.setValue("isotope:stdev", 1.0, "Standard deviation of gaussian applied to the averagine isotopic pattern to simulate the inaccuracy of the mass spectrometer.", {"advanced"});
    defaults_.setValue("isotope:maximum", 100, "Maximum isotopic rank to be considered.", {"advanced"});
    defaults_.setValue("interpolation_step", 0.1, "Sampling rate for the interpolation of the model function.", {"advanced"});

    defaultsToParam_();
  }

  IsotopeFitter1D::IsotopeFitter1D(const IsotopeFitter1D& source) :
    MaxLikeliFitter1D(source)
  {
    updateMembers_();
  }

  IsotopeFitter1D::~IsotopeFitter1D() = default;

  IsotopeFitter1D& IsotopeFitter1D::operator=(const IsotopeFitter1D& source)
  {
    if (&source == this)
    {
      return *this;
    }
    MaxLikeliFitter1D::operator=(source);
    updateMembers_();

    return *this;
  }

  IsotopeFitter1D::QualityType IsotopeFitter1D::fit1d(const RawDataArrayType& set, std::unique_ptr<InterpolationModel>& model)
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
    if (charge_ == 0)
    {
  model = std::unique_ptr<InterpolationModel>(new GaussModel());
      model->setInterpolationStep(interpolation_step_);

      Param tmp;
      tmp.setValue("bounding_box:min", min_bb);
      tmp.setValue("bounding_box:max", max_bb);
      tmp.setValue("statistics:variance", statistics_.variance());
      tmp.setValue("statistics:mean", statistics_.mean());
      model->setParameters(tmp);
    }
    else
    {
      model = std::unique_ptr<InterpolationModel>(new IsotopeModel());

      Param iso_param = this->param_.copy("isotope_model:", true);
      iso_param.removeAll("stdev");
      model->setParameters(iso_param);
      model->setInterpolationStep(interpolation_step_);

      Param tmp;
      tmp.setValue("statistics:mean", statistics_.mean());
      tmp.setValue("charge", static_cast<Int>(charge_));
      tmp.setValue("isotope:mode:GaussianSD", isotope_stdev_);
      tmp.setValue("isotope:maximum", max_isotope_);

      model->setParameters(tmp);
      (dynamic_cast<IsotopeModel*>(model.get()))->setSamples((dynamic_cast<IsotopeModel*>(model.get()))->getFormula());
    }

    // fit offset
    QualityType quality;
    quality = fitOffset_(model, set, stdev, stdev, interpolation_step_);
    if (std::isnan(quality))
    {
      quality = -1.0;
    }
    return quality;
  }

  void IsotopeFitter1D::updateMembers_()
  {
    MaxLikeliFitter1D::updateMembers_();
    statistics_.setVariance(param_.getValue("statistics:variance"));
    charge_ = param_.getValue("charge");
    isotope_stdev_ = param_.getValue("isotope:stdev");
    max_isotope_ = param_.getValue("isotope:maximum");
  }

}
