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

  void PeakIntegrator::estimateBackground(const MSChromatogram& chromatogram, const double& left, const double& right)
  {
    estimateBackground_(chromatogram, left, right);
  }

  void PeakIntegrator::estimateBackground(const MSSpectrum& spectrum, const double& left, const double& right)
  {
    estimateBackground_(spectrum, left, right);
  }

  void PeakIntegrator::integratePeak(const MSChromatogram& chromatogram, const double& left, const double& right)
  {
    integratePeak_(chromatogram, left, right);
  }

  void PeakIntegrator::integratePeak(const MSSpectrum& spectrum, const double& left, const double& right)
  {
    integratePeak_(spectrum, left, right);
  }

  double PeakIntegrator::simpson(MSChromatogram::ConstIterator it_begin, MSChromatogram::ConstIterator it_end) const
  {
    return simpson_(it_begin, it_end);
  }

  double PeakIntegrator::simpson(MSSpectrum::ConstIterator it_begin, MSSpectrum::ConstIterator it_end) const
  {
    return simpson_(it_begin, it_end);
  }

  void PeakIntegrator::calculatePeakShapeMetrics(const MSChromatogram& chromatogram, const double& left, const double& right)
  {
    calculatePeakShapeMetrics_(chromatogram, left, right);
  }

  void PeakIntegrator::calculatePeakShapeMetrics(const MSSpectrum& spectrum, const double& left, const double& right)
  {
    calculatePeakShapeMetrics_(spectrum, left, right);
  }

  double PeakIntegrator::getPeakArea() const { return peak_area_; }
  double PeakIntegrator::getPeakHeight() const { return peak_height_; }
  double PeakIntegrator::getPeakApexPos() const { return peak_apex_pos_; }
  double PeakIntegrator::getBackgroundHeight() const { return background_height_; }
  double PeakIntegrator::getBackgroundArea() const { return background_area_; }
  PeakIntegrator::PeakShapeMetrics PeakIntegrator::getPeakShapeMetrics() const { return psm_; }

  void PeakIntegrator::getDefaultParameters(Param& params)
  {
    params.clear();
    params.setValue("integration_type", "intensity_sum", "The integration technique to use in integratePeak() and estimateBackground().");
    params.setValidStrings("integration_type", ListUtils::create<String>("intensity_sum,simpson,trapezoid"));
    params.setValue("baseline_type", "base_to_base", "The baseline type to use in estimateBackground().");
    params.setValidStrings("baseline_type", ListUtils::create<String>("base_to_base,vertical_division"));
  }

  void PeakIntegrator::updateMembers_()
  {
    integration_type_ = (String)param_.getValue("integration_type");
    baseline_type_ = (String)param_.getValue("baseline_type");
  }

  template <typename PeakContainerT>
  void PeakIntegrator::estimateBackground_(const PeakContainerT& p, const double& left, const double& right)
  {
    const double int_l = p.PosBegin(left)->getIntensity();
    const double int_r = (p.PosEnd(right)-1)->getIntensity();
    const double delta_int = int_r - int_l;
    const double delta_pos = (p.PosEnd(right)-1)->getPos() - p.PosBegin(left)->getPos();
    const double min_int_pos = int_r <= int_l ? (p.PosEnd(right)-1)->getPos() : p.PosBegin(left)->getPos();
    const double delta_int_apex = std::fabs(delta_int) * std::fabs(min_int_pos - peak_apex_pos_) / delta_pos;
    background_height_ = std::min(int_r, int_l) + delta_int_apex;
    double background = 0.0;
    if (baseline_type_ == "base_to_base")
    {
      if (integration_type_ == "trapezoid" || integration_type_ == "simpson")
      {
        // formula for calculating the background using the trapezoidal rule
        // background = intensity_min*delta_pos + 0.5*delta_int*delta_pos;
        background = delta_pos * (std::min(int_r, int_l) + 0.5 * std::fabs(delta_int));
      }
      else if (integration_type_ == "intensity_sum")
      {
        // calculate the background using the formula
        // y = mx + b where x = rt or mz, m = slope, b = left intensity
        // sign of delta_int will determine line direction
        // background += delta_int / delta_pos * (it->getPos() - left) + int_l;
        UInt n_points = 0;
        for (auto it=p.PosBegin(left); it!=p.PosEnd(right); ++it, ++n_points)
        {
          background += it->getPos();
        }
        background = (background - n_points * p.PosBegin(left)->getPos()) * delta_int / delta_pos + n_points * int_l;
      }
    }
    else if (baseline_type_ == "vertical_division")
    {
      if (integration_type_ == "trapezoid" || integration_type_ == "simpson")
      {
        background = delta_pos * std::min(int_r, int_l);
      }
      else if (integration_type_ == "intensity_sum")
      {
        UInt n_points = 0;
        for (auto it=p.PosBegin(left); it!=p.PosEnd(right); ++it, ++n_points)
        { }
        background = std::min(int_r, int_l) * n_points;
      }
    }
    background_area_ = background;
  }

  template <typename PeakContainerT>
  void PeakIntegrator::integratePeak_(PeakContainerT p, const double& left, const double& right)
  {
    peak_area_ = 0.0;
    peak_height_ = -1.0;
    peak_apex_pos_ = -1.0;
    UInt n_points = 0;
    for (auto it=p.PosBegin(left); it!=p.PosEnd(right); ++it, ++n_points)
    {
      if (peak_height_ < it->getIntensity())
      {
        peak_height_ = it->getIntensity();
        peak_apex_pos_ = it->getPos();
      }
    }

    if (integration_type_ == "trapezoid")
    {
      for (auto it=p.PosBegin(left); it!=p.PosEnd(right)-1; ++it)
      {
        peak_area_ += ((it+1)->getPos() - it->getPos()) * ((it->getIntensity() + (it+1)->getIntensity()) / 2.0);
      }
    }
    else if (integration_type_ == "simpson")
    {
      if (n_points < 3)
      {
        LOG_DEBUG << std::endl << "Error in integratePeak: number of points must be >=3 for Simpson's rule" << std::endl;
        return;
      }
      if (n_points % 2)
      {
        peak_area_ = simpson(p.PosBegin(left), p.PosEnd(right));
      }
      else
      {
        double areas[4] = {-1.0, -1.0, -1.0, -1.0};
        areas[0] = simpson(p.PosBegin(left), p.PosEnd(right) - 1);   // without last point
        areas[1] = simpson(p.PosBegin(left) + 1, p.PosEnd(right));   // without first point
        if (p.begin() <= p.PosBegin(left) - 1)
        {
          areas[2] = simpson(p.PosBegin(left) - 1, p.PosEnd(right)); // with one more point on the left
        }
        if (p.PosEnd(right) < p.end())
        {
          areas[3] = simpson(p.PosBegin(left), p.PosEnd(right) + 1); // with one more point on the right
        }
        UInt valids = 0;
        for (auto area : areas)
        {
          if (area != -1.0)
          {
            peak_area_ += area;
            ++valids;
          }
        }
        peak_area_ /= valids;
      }
    }
    else if (integration_type_ == "intensity_sum")
    {
      std::cout << "\nWARNING: intensity_sum method is being used.\n";
      for (auto it=p.PosBegin(left); it!=p.PosEnd(right); ++it)
      {
        peak_area_ += it->getIntensity();
      }
    }
    else
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Please set a valid value for the parameter \"integration_type\".\n");
    }
  }

  template <typename PeakContainerConstIteratorT>
  double PeakIntegrator::simpson_(PeakContainerConstIteratorT it_begin, PeakContainerConstIteratorT it_end) const
  {
    double integral = 0.0;
    for (auto it=it_begin+1; it<it_end-1; it=it+2)
    {
      const double h = it->getPos() - (it-1)->getPos();
      const double k = (it+1)->getPos() - it->getPos();
      const double y_h = (it-1)->getIntensity();
      const double y_0 = it->getIntensity();
      const double y_k = (it+1)->getIntensity();
      integral += (1.0/6.0) * (h+k) * ((2.0-k/h)*y_h + (pow(h+k,2)/(h*k))*y_0 + (2.0-h/k)*y_k);
    }
    return integral;
  }

  template <typename PeakContainerT>
  void PeakIntegrator::calculatePeakShapeMetrics_(const PeakContainerT& p, const double& left, const double& right)
  {
    psm_.points_across_baseline = 0;
    psm_.points_across_half_height = 0;
    for (auto it = p.PosBegin(left); it != p.PosEnd(right); ++it)
    {
      const double intensity = it->getIntensity();
      const double intensity_prev = (it-1)->getIntensity();
      const double position = it->getPos();
      const double position_prev = (it-1)->getPos();
      // start and end positions (rt or mz)
      if (position < peak_apex_pos_ && psm_.points_across_baseline > 1) // start positions
      {
        const double d_int_times_d_pos = (intensity - intensity_prev) * (position - position_prev);
        if (intensity >= 0.05 * peak_height_ && intensity_prev < 0.05 * peak_height_)
        {
          const double height_5 = intensity - 0.05 * peak_height_;
          psm_.start_position_at_5 = position - d_int_times_d_pos / height_5;
        }
        if (intensity >= 0.1 * peak_height_ && intensity_prev < 0.1 * peak_height_)
        {
          const double height_10 = intensity - 0.1 * peak_height_;
          psm_.start_position_at_10 = position - d_int_times_d_pos / height_10;
        }
        if (intensity >= 0.5 * peak_height_ && intensity_prev < 0.5 * peak_height_)
        {
          const double height_50 = intensity - 0.5 * peak_height_;
          psm_.start_position_at_50 = position - d_int_times_d_pos / height_50;
        }
      }
      else if (position > peak_apex_pos_) // end positions
      {
        const double d_int_times_d_pos = (intensity_prev - intensity) * (position - position_prev);
        if (intensity <= 0.05 * peak_height_ && intensity_prev > 0.05 * peak_height_)
        {
          const double height_5 = 0.05 * peak_height_ - intensity;
          psm_.end_position_at_5 = position - d_int_times_d_pos / height_5;
        }
        if (intensity <= 0.1 * peak_height_ && intensity_prev > 0.1 * peak_height_)
        {
          const double height_10 = 0.1 * peak_height_ - intensity;
          psm_.end_position_at_10 = position - d_int_times_d_pos / height_10;
        }
        if (intensity <= 0.5 * peak_height_ && intensity_prev > 0.5 * peak_height_)
        {
          const double height_50 = 0.5 * peak_height_ - intensity;
          psm_.end_position_at_50 = position - d_int_times_d_pos / height_50;
        }
      }
      // points across the peak
      ++(psm_.points_across_baseline);
      if (intensity >= 0.5 * peak_height_)
      {
        ++(psm_.points_across_half_height);
      }
    }

    // peak widths
    psm_.width_at_5 = psm_.end_position_at_5 - psm_.start_position_at_5;
    psm_.width_at_10 = psm_.end_position_at_10 - psm_.start_position_at_10;
    psm_.width_at_50 = psm_.end_position_at_50 - psm_.start_position_at_50;
    psm_.total_width = (p.PosEnd(right)-1)->getPos() - p.PosBegin(left)->getPos();
    psm_.slope_of_baseline = (p.PosEnd(right)-1)->getIntensity() - p.PosBegin(left)->getIntensity();
    psm_.baseline_delta_2_height = psm_.slope_of_baseline / peak_height_;

    // other
    psm_.tailing_factor = psm_.width_at_5 / std::min(peak_apex_pos_ - psm_.start_position_at_5, psm_.end_position_at_5 - peak_apex_pos_);
    psm_.asymmetry_factor = std::min(peak_apex_pos_ - psm_.start_position_at_10, psm_.end_position_at_10 - peak_apex_pos_) /
      std::max(peak_apex_pos_ - psm_.start_position_at_10, psm_.end_position_at_10 - peak_apex_pos_);
  }
}
