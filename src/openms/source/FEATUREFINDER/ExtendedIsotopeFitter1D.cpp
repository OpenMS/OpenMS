// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/ExtendedIsotopeFitter1D.h>

#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/FEATUREFINDER/ExtendedIsotopeModel.h>

namespace OpenMS
{

  ExtendedIsotopeFitter1D::ExtendedIsotopeFitter1D() :
    MaxLikeliFitter1D()
  {
    setName("ExtendedIsotopeFitter1D");

    defaults_.setValue("statistics:variance", 1.0, "Variance of the model.", {"advanced"});
    defaults_.setValue("charge", 1, "Charge state of the model.", {"advanced"});
    defaults_.setValue("isotope:stdev", 0.1, "Standard deviation of gaussian applied to the averagine isotopic pattern to simulate the inaccuracy of the mass spectrometer.", {"advanced"});
    defaults_.setValue("isotope:monoisotopic_mz", 1.0, "Monoisotopic m/z of the model.", {"advanced"});
    defaults_.setValue("isotope:maximum", 100, "Maximum isotopic rank to be considered.", {"advanced"});
    defaults_.setValue("interpolation_step", 0.2, "Sampling rate for the interpolation of the model function.", {"advanced"});

    defaultsToParam_();
  }

  ExtendedIsotopeFitter1D::ExtendedIsotopeFitter1D(const ExtendedIsotopeFitter1D& source) :
    MaxLikeliFitter1D(source)
  {
    updateMembers_();
  }

  ExtendedIsotopeFitter1D::~ExtendedIsotopeFitter1D() = default;

  ExtendedIsotopeFitter1D& ExtendedIsotopeFitter1D::operator=(const ExtendedIsotopeFitter1D& source)
  {
    if (&source == this)
    {
      return *this;
    }
    MaxLikeliFitter1D::operator=(source);
    updateMembers_();

    return *this;
  }

  ExtendedIsotopeFitter1D::QualityType ExtendedIsotopeFitter1D::fit1d(const RawDataArrayType& set, std::unique_ptr<InterpolationModel>& model)
  {
    // build model
    if (charge_ == 0)
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
      model = std::unique_ptr<InterpolationModel>(new ExtendedIsotopeModel());

      Param iso_param = this->param_.copy("isotope_model:", true);
      iso_param.removeAll("stdev");
      model->setParameters(iso_param);
      model->setInterpolationStep(interpolation_step_);

      Param tmp;
      tmp.setValue("isotope:monoisotopic_mz", monoisotopic_mz_);
      tmp.setValue("charge", static_cast<Int>(charge_));
      tmp.setValue("isotope:stdev", isotope_stdev_);
      tmp.setValue("isotope:maximum", max_isotope_);
      model->setParameters(tmp);
    }


    // calculate pearson correlation
    std::vector<float> real_data;
    real_data.reserve(set.size());
    std::vector<float> model_data;
    model_data.reserve(set.size());

    for (Size i = 0; i < set.size(); ++i)
    {
      real_data.push_back(set[i].getIntensity());
      model_data.push_back(model->getIntensity(DPosition<1>(set[i].getPosition())));
    }

    QualityType correlation = Math::pearsonCorrelationCoefficient(real_data.begin(), real_data.end(), model_data.begin(), model_data.end());
    if (std::isnan(correlation))
    {
      correlation = -1.0;
    }
    return correlation;
  }

  void ExtendedIsotopeFitter1D::updateMembers_()
  {
    MaxLikeliFitter1D::updateMembers_();
    statistics_.setVariance(param_.getValue("statistics:variance"));
    charge_ = param_.getValue("charge");
    isotope_stdev_ = param_.getValue("isotope:stdev");
    monoisotopic_mz_ = param_.getValue("isotope:monoisotopic_mz");
    max_isotope_ = param_.getValue("isotope:maximum");

  }

}
