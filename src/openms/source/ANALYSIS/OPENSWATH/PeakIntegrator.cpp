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
  double PeakIntegrator::getWidthAt5() const { return width_at_5_; }
  double PeakIntegrator::getWidthAt10() const { return width_at_10_; }
  double PeakIntegrator::getWidthAt50() const { return width_at_50_; }
  double PeakIntegrator::getStartTimeAt5() const { return start_time_at_5_; }
  double PeakIntegrator::getStartTimeAt10() const { return start_time_at_10_; }
  double PeakIntegrator::getStartTimeAt50() const { return start_time_at_50_; }
  double PeakIntegrator::getEndTimeAt5() const { return end_time_at_5_; }
  double PeakIntegrator::getEndTimeAt10() const { return end_time_at_10_; }
  double PeakIntegrator::getEndTimeAt50() const { return end_time_at_50_; }
  double PeakIntegrator::getTotalWidth() const { return total_width_; }
  double PeakIntegrator::getTailingFactor() const { return tailing_factor_; }
  double PeakIntegrator::getAsymmetryFactor() const { return asymmetry_factor_; }
  double PeakIntegrator::getBaselineDeltaToHeight() const { return baseline_delta_2_height_; }
  double PeakIntegrator::getSlopeOfBaseline() const { return slope_of_baseline_; }
  Int PeakIntegrator::getPointsAcrossBaseline() const { return points_across_baseline_; }
  Int PeakIntegrator::getPointsAcrossHalfHeight() const { return points_across_half_height_; }

  const std::map<String, double> PeakIntegrator::getPeakShapeMetrics() const
  {
    std::map<String, double> m
    {
      { "width_at_5", width_at_5_ },
      { "width_at_10", width_at_10_ },
      { "width_at_50", width_at_50_ },
      { "start_time_at_5", start_time_at_5_ },
      { "start_time_at_10", start_time_at_10_ },
      { "start_time_at_50", start_time_at_50_ },
      { "end_time_at_5", end_time_at_5_ },
      { "end_time_at_10", end_time_at_10_ },
      { "end_time_at_50", end_time_at_50_ },
      { "total_width", total_width_ },
      { "tailing_factor", tailing_factor_ },
      { "asymmetry_factor", asymmetry_factor_ },
      { "baseline_delta_to_height", baseline_delta_2_height_ },
      { "slope_of_baseline", slope_of_baseline_ },
      { "points_across_baseline", points_across_baseline_ },
      { "points_across_half_height", points_across_half_height_ },
    };
    return m;
  }

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
          ;
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
    else
    {
      std::cout << std::endl << "WARNING: intensity_sum method is being used." << std::endl;
      for (auto it=p.PosBegin(left); it!=p.PosEnd(right); ++it)
      {
        peak_area_ += it->getIntensity();
      }
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
    points_across_baseline_ = 0;
    points_across_half_height_ = 0;
    for (auto it = p.PosBegin(left); it != p.PosEnd(right); ++it)
    {
      const double intensity = it->getIntensity();
      const double intensity_prev = (it-1)->getIntensity();
      const double position = it->getPos();
      const double position_prev = (it-1)->getPos();
      // start and end positions (rt or mz)
      if (position < peak_apex_pos_ && points_across_baseline_ > 1) // start positions
      {
        const double d_int_times_d_rt = (intensity - intensity_prev) * (position - position_prev);
        if (intensity >= 0.05 * peak_height_ && intensity_prev < 0.05 * peak_height_)
        {
          const double height_5 = intensity - 0.05 * peak_height_;
          start_time_at_5_ = position - d_int_times_d_rt / height_5;
        }
        if (intensity >= 0.1 * peak_height_ && intensity_prev < 0.1 * peak_height_)
        {
          const double height_10 = intensity - 0.1 * peak_height_;
          start_time_at_10_ = position - d_int_times_d_rt / height_10;
        }
        if (intensity >= 0.5 * peak_height_ && intensity_prev < 0.5 * peak_height_)
        {
          const double height_50 = intensity - 0.5 * peak_height_;
          start_time_at_50_ = position - d_int_times_d_rt / height_50;
        }
      }
      else if (position > peak_apex_pos_) // end positions
      {
        const double d_int_times_d_rt = (intensity_prev - intensity) * (position - position_prev);
        if (intensity <= 0.05 * peak_height_ && intensity_prev > 0.05 * peak_height_)
        {
          const double height_5 = 0.05 * peak_height_ - intensity;
          end_time_at_5_ = position - d_int_times_d_rt / height_5;
        }
        if (intensity <= 0.1 * peak_height_ && intensity_prev > 0.1 * peak_height_)
        {
          const double height_10 = 0.1 * peak_height_ - intensity;
          end_time_at_10_ = position - d_int_times_d_rt / height_10;
        }
        if (intensity <= 0.5 * peak_height_ && intensity_prev > 0.5 * peak_height_)
        {
          const double height_50 = 0.5 * peak_height_ - intensity;
          end_time_at_50_ = position - d_int_times_d_rt / height_50;
        }
      }
      // points across the peak
      ++points_across_baseline_;
      if (intensity >= 0.5 * peak_height_)
      {
        ++points_across_half_height_;
      }
    }

    // peak widths
    width_at_5_ = end_time_at_5_ - start_time_at_5_;
    width_at_10_ = end_time_at_10_ - start_time_at_10_;
    width_at_50_ = end_time_at_50_ - start_time_at_50_;
    total_width_ = (p.PosEnd(right)-1)->getPos() - p.PosBegin(left)->getPos();
    slope_of_baseline_ = (p.PosEnd(right)-1)->getIntensity() - p.PosBegin(left)->getIntensity();
    baseline_delta_2_height_ = slope_of_baseline_ / peak_height_;

    // other
    tailing_factor_ = width_at_5_ / std::min(peak_apex_pos_ - start_time_at_5_, end_time_at_5_ - peak_apex_pos_);
    asymmetry_factor_ = std::min(peak_apex_pos_ - start_time_at_10_, end_time_at_10_ - peak_apex_pos_) /
      std::max(peak_apex_pos_ - start_time_at_10_, end_time_at_10_ - peak_apex_pos_);
  }
}
