// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeFitter1D.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/Factory.h>

#include <numeric>
#include <boost/math/special_functions/fpclassify.hpp>

namespace OpenMS
{

  ExtendedIsotopeFitter1D::ExtendedIsotopeFitter1D() :
    MaxLikeliFitter1D()
  {
    setName(getProductName());

    defaults_.setValue("statistics:variance", 1.0, "Variance of the model.", ListUtils::create<String>("advanced"));
    defaults_.setValue("charge", 1, "Charge state of the model.", ListUtils::create<String>("advanced"));
    defaults_.setValue("isotope:stdev", 0.1, "Standard deviation of gaussian applied to the averagine isotopic pattern to simulate the inaccuracy of the mass spectrometer.", ListUtils::create<String>("advanced"));
    defaults_.setValue("isotope:monoisotopic_mz", 1.0, "Monoisotopic m/z of the model.", ListUtils::create<String>("advanced"));
    defaults_.setValue("isotope:maximum", 100, "Maximum isotopic rank to be considered.", ListUtils::create<String>("advanced"));
    defaults_.setValue("interpolation_step", 0.2, "Sampling rate for the interpolation of the model function.", ListUtils::create<String>("advanced"));

    defaultsToParam_();
  }

  ExtendedIsotopeFitter1D::ExtendedIsotopeFitter1D(const ExtendedIsotopeFitter1D & source) :
    MaxLikeliFitter1D(source)
  {
    updateMembers_();
  }

  ExtendedIsotopeFitter1D::~ExtendedIsotopeFitter1D()
  {
  }

  ExtendedIsotopeFitter1D & ExtendedIsotopeFitter1D::operator=(const ExtendedIsotopeFitter1D & source)
  {
    if (&source == this)
      return *this;

    MaxLikeliFitter1D::operator=(source);
    updateMembers_();

    return *this;
  }

  ExtendedIsotopeFitter1D::QualityType ExtendedIsotopeFitter1D::fit1d(const RawDataArrayType & set, InterpolationModel * & model)
  {
    // build model
    if (charge_ == 0)
    {
      // Calculate bounding box
      min_ = max_ = set[0].getPos();
      for (UInt pos = 1; pos < set.size(); ++pos)
      {
        CoordinateType tmp = set[pos].getPos();
        if (min_ > tmp)
          min_ = tmp;
        if (max_ < tmp)
          max_ = tmp;
      }

      // Enlarge the bounding box by a few multiples of the standard deviation
      {
        stdev1_ = sqrt(statistics_.variance()) * tolerance_stdev_box_;
        min_ -= stdev1_;
        max_ += stdev1_;
      }

      model = static_cast<InterpolationModel *>(Factory<BaseModel<1> >::create("GaussModel"));
      model->setInterpolationStep(interpolation_step_);

      Param tmp;
      tmp.setValue("bounding_box:min", min_);
      tmp.setValue("bounding_box:max", max_);
      tmp.setValue("statistics:variance", statistics_.variance());
      tmp.setValue("statistics:mean", statistics_.mean());
      model->setParameters(tmp);
    }
    else
    {
      model = static_cast<InterpolationModel *>(Factory<BaseModel<1> >::create("ExtendedIsotopeModel"));

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
    if (boost::math::isnan(correlation))
      correlation = -1.0;

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
