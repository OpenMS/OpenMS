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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>

namespace OpenMS
{
  PeakIntegrator::PeakIntegrator() :
    DefaultParamHandler("PeakIntegrator")
  {
    getDefaultParameters(defaults_);
    defaultsToParam_(); // write defaults into Param object param_
  }

  PeakIntegrator::~PeakIntegrator() {}

  double PeakIntegrator::estimateBackground(
    const MSChromatogram& chromatogram,
    const double& left,
    const double& right
  )
  {
    const double intensity_l = chromatogram.RTBegin(left)->getIntensity();
    const double intensity_r = (chromatogram.RTEnd(right)-1)->getIntensity();
    const double delta_int = intensity_r - intensity_l;
    const double delta_rt = (chromatogram.RTEnd(right)-1)->getRT() - chromatogram.RTBegin(left)->getRT();
    double background = 0.0;
    if (integration_type_ == "trapezoid" || integration_type_ == "simpson")
    {
      // formula for calculating the background using the trapezoidal rule
      // background = intensity_min*delta_rt + 0.5*delta_int*delta_rt;
      background = delta_rt * (std::min(intensity_r, intensity_l) + 0.5 * std::fabs(delta_int));
    }
    else
    {
      // calculate the background using the formula
      // y = mx + b where x = retention time, m = slope, b = left intensity
      // sign of delta_int will determine line direction
      // background += delta_int / delta_rt * (it->getRT() - left) + intensity_l;
      UInt i = 0;
      for (auto it=chromatogram.RTBegin(left); it!=chromatogram.RTEnd(right); ++it, ++i)
      {
        background += it->getRT();
      }
      background = (background - i * chromatogram.RTBegin(left)->getRT()) * delta_int / delta_rt + i * intensity_l;
    }
    return background;
  }

  void PeakIntegrator::integratePeak(
    const MSChromatogram& chromatogram,
    const double& left,
    const double& right,
    double& peak_area,
    double& peak_height,
    double& peak_apex_pos
  )
  {
    peak_area = 0.0;
    peak_height = -1.0;
    peak_apex_pos = -1.0;
    UInt n_points = 0;
    for (auto it=chromatogram.RTBegin(left); it!=chromatogram.RTEnd(right); ++it, ++n_points)
      ;

    if (getIntegrationType() == "trapezoid")
    {
      for (auto it=chromatogram.RTBegin(left); it!=chromatogram.RTEnd(right)-1; ++it)
      {
        peak_area += ((it+1)->getRT() - it->getRT()) * ((it->getIntensity() + (it+1)->getIntensity()) / 2.0);
        if (peak_height < it->getIntensity())
        {
          peak_height = it->getIntensity();
          peak_apex_pos = it->getRT();
        }
      }
    }
    else if (getIntegrationType() == "simpson")
    {
      if (n_points < 3 || !(n_points % 2))
      {
        LOG_DEBUG << std::endl << "Error in integratePeak: number of points must be >=3 and odd for Simpson's rule" << std::endl;
        return;
      }
      for (auto it=chromatogram.RTBegin(left); it<chromatogram.RTEnd(right)-2; it=it+2)
      {
        double h = (it+1)->getRT() - it->getRT();
        double k = (it+2)->getRT() - (it+1)->getRT();
        double y_h = it->getIntensity();
        double y_0 = (it+1)->getIntensity();
        double y_k = (it+2)->getIntensity();
        peak_area += (1.0/6.0) * (h+k) * ((2.0-k/h)*y_h + (pow(h+k,2)/(h*k))*y_0 + (2.0-h/k)*y_k);
        if (peak_height < it->getIntensity())
        {
          peak_height = it->getIntensity();
          peak_apex_pos = it->getRT();
        }
      }
    }
    else
    {
      std::cout << std::endl << "WARNING: intensity_sum method is being used." << std::endl;
      double intensity_sum = 0.0;
      for (auto it=chromatogram.RTBegin(left); it!=chromatogram.RTEnd(right); ++it)
      {
        intensity_sum += it->getIntensity();
        if (peak_height < it->getIntensity())
        {
          peak_height = it->getIntensity();
          peak_apex_pos = it->getRT();
        }
      }
      peak_area = intensity_sum / n_points;
    }
  }

  void PeakIntegrator::setIntegrationType(const String& integration_type)
  {
    integration_type_ = integration_type;
  }

  String PeakIntegrator::getIntegrationType() const
  {
    return integration_type_;
  }

  void PeakIntegrator::setBaselineType(const String& baseline_type)
  {
    baseline_type_ = baseline_type;
  }

  String PeakIntegrator::getBaselineType() const
  {
    return baseline_type_;
  }

  void PeakIntegrator::setPeakModel(const String& peak_model)
  {
    peak_model_ = peak_model;
  }

  String PeakIntegrator::getPeakModel() const
  {
    return peak_model_;
  }

  void PeakIntegrator::getDefaultParameters(Param& params)
  {
    params.clear();
    // TODO improve descriptions
    params.setValue("integration_type", "trapezoid", "Integration method to use.");
    params.setValidStrings("integration_type", ListUtils::create<String>("intensity_sum,simpson,trapezoid"));

    params.setValue("baseline_type", "vertical_division", "Type of baseline to use.");
    params.setValidStrings("baseline_type", ListUtils::create<String>("base_to_base,vertical_division"));

    params.setValue("peak_model", "none", "Peak model.");
    params.setValidStrings("peak_model", ListUtils::create<String>("none"));
  }

  void PeakIntegrator::updateMembers_()
  {
    integration_type_ = (String)param_.getValue("integration_type");
    baseline_type_ = (String)param_.getValue("baseline_type");
    peak_model_ = (String)param_.getValue("peak_model");
  }
}
