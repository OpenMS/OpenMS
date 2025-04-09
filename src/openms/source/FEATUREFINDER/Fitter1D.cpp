// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/Fitter1D.h>

// include derived classes here
#include <OpenMS/FEATUREFINDER/GaussFitter1D.h>
#include <OpenMS/FEATUREFINDER/BiGaussFitter1D.h>
#include <OpenMS/FEATUREFINDER/IsotopeFitter1D.h>
#include <OpenMS/FEATUREFINDER/ExtendedIsotopeFitter1D.h>
#include <OpenMS/FEATUREFINDER/EmgFitter1D.h>

namespace OpenMS
{
  Fitter1D::Fitter1D() :
    DefaultParamHandler("Fitter1D")
  {
    defaults_.setValue("interpolation_step", 0.2, "Sampling rate for the interpolation of the model function.", {"advanced"});
    defaults_.setValue("statistics:mean", 1.0, "Centroid position of the model.", {"advanced"});
    defaults_.setValue("statistics:variance", 1.0, "The variance of the model.", {"advanced"});
    defaults_.setValue("tolerance_stdev_bounding_box", 3.0, "Bounding box has range [minimim of data, maximum of data] enlarged by tolerance_stdev_bounding_box times the standard deviation of the data.", {"advanced"});

    defaultsToParam_();
  }

  Fitter1D::Fitter1D(const Fitter1D& source) :
    DefaultParamHandler(source)
  {
    setParameters(source.getParameters());
    updateMembers_();
  }

  Fitter1D::~Fitter1D() = default;

  Fitter1D& Fitter1D::operator=(const Fitter1D& source)
  {
    if (&source == this)
    {
      return *this;
    }
    DefaultParamHandler::operator=(source);
    setParameters(source.getParameters());
    updateMembers_();

    return *this;
  }

  void Fitter1D::updateMembers_()
  {
    tolerance_stdev_box_ = param_.getValue("tolerance_stdev_bounding_box");
    interpolation_step_ = param_.getValue("interpolation_step");
    statistics_.setMean(param_.getValue("statistics:mean"));
    statistics_.setVariance(param_.getValue("statistics:variance"));
  }

  Fitter1D::QualityType Fitter1D::fit1d(const RawDataArrayType& /* range */, std::unique_ptr<InterpolationModel>& /* model */)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

} // namespace OpenMS
