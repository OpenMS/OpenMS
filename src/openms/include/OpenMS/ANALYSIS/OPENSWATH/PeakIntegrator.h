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

#ifndef OPENMS_ANALYSIS_OPENSWATH_PEAKINTEGRATOR_H
#define OPENMS_ANALYSIS_OPENSWATH_PEAKINTEGRATOR_H

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  class OPENMS_DLLAPI PeakIntegrator :
    public DefaultParamHandler
  {
public:
    PeakIntegrator();
    virtual ~PeakIntegrator();

    void estimateBackground(const MSChromatogram& chromatogram, const double& left, const double& right);
    void estimateBackground(const MSSpectrum& spectrum, const double& left, const double& right);

    void integratePeak(const MSChromatogram& chromatogram, const double& left, const double& right);
    void integratePeak(const MSSpectrum& spectrum, const double& left, const double& right);

    void calculatePeakShapeMetrics(const MSChromatogram& chromatogram, const double& left, const double& right);
    void calculatePeakShapeMetrics(const MSSpectrum& spectrum, const double& left, const double& right);

    double getPeakArea() const;
    double getPeakHeight() const;
    double getPeakApexRT() const;
    double getBackgroundHeight() const;
    double getBackgroundArea() const;
    double getWidthAt5() const;
    double getWidthAt10() const;
    double getWidthAt50() const;
    double getStartTimeAt5() const;
    double getStartTimeAt10() const;
    double getStartTimeAt50() const;
    double getEndTimeAt5() const;
    double getEndTimeAt10() const;
    double getEndTimeAt50() const;
    double getTotalWidth() const;
    double getTailingFactor() const;
    double getAsymmetryFactor() const;
    double getBaselineDeltaToHeight() const;
    double getSlopeOfBaseline() const;
    Int getPointsAcrossBaseline() const;
    Int getPointsAcrossHalfHeight() const;

    const std::map<String, double> getPeakShapeMetrics() const;

    void getDefaultParameters(Param& params);

protected:
    void updateMembers_();

    template<class PeakContainerT>
    void estimateBackground_(const PeakContainerT& p, const double& left, const double& right)
    {
      const double int_l = p.PosBegin(left)->getIntensity();
      const double int_r = (p.PosEnd(right)-1)->getIntensity();
      const double delta_int = int_r - int_l;
      const double delta_rt = (p.PosEnd(right)-1)->getPos() - p.PosBegin(left)->getPos();
      const double rt_min = int_r <= int_l ? (p.PosEnd(right)-1)->getPos() : p.PosBegin(left)->getPos();
      const double delta_int_apex = std::fabs(delta_int) * std::fabs(rt_min - peak_apex_rt_) / delta_rt;
      background_height_ = std::min(int_r, int_l) + delta_int_apex;
      double background = 0.0;
      if (baseline_type_ == "base_to_base")
      {
        if (integration_type_ == "trapezoid" || integration_type_ == "simpson")
        {
          // formula for calculating the background using the trapezoidal rule
          // background = intensity_min*delta_rt + 0.5*delta_int*delta_rt;
          background = delta_rt * (std::min(int_r, int_l) + 0.5 * std::fabs(delta_int));
        }
        else if (integration_type_ == "intensity_sum")
        {
          // calculate the background using the formula
          // y = mx + b where x = retention time, m = slope, b = left intensity
          // sign of delta_int will determine line direction
          // background += delta_int / delta_rt * (it->getPos() - left) + int_l;
          UInt n_points = 0;
          for (auto it=p.PosBegin(left); it!=p.PosEnd(right); ++it, ++n_points)
          {
            background += it->getPos();
          }
          background = (background - n_points * p.PosBegin(left)->getPos()) * delta_int / delta_rt + n_points * int_l;
        }
      }
      else if (baseline_type_ == "vertical_division")
      {
        if (integration_type_ == "trapezoid" || integration_type_ == "simpson")
        {
          background = delta_rt * std::min(int_r, int_l);
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

    template<class PeakContainerT>
    void integratePeak_(PeakContainerT p, const double& left, const double& right)
    {
      peak_area_ = 0.0;
      peak_height_ = -1.0;
      peak_apex_rt_ = -1.0;
      UInt n_points = 0;
      for (auto it=p.PosBegin(left); it!=p.PosEnd(right); ++it, ++n_points)
      {
        if (peak_height_ < it->getIntensity())
        {
          peak_height_ = it->getIntensity();
          peak_apex_rt_ = it->getPos();
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
          double areas[4] = {};
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
            if (area)
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

    template<class PeakContainerConstIteratorT>
    double simpson_(PeakContainerConstIteratorT it_begin, PeakContainerConstIteratorT it_end) const
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

    template<class PeakContainerT>
    void calculatePeakShapeMetrics_(const PeakContainerT& p, const double& left, const double& right)
    {
      points_across_baseline_ = 0;
      points_across_half_height_ = 0;
      double start_intensity(0), end_intensity(0);
      double delta_rt, delta_int, height_5, height_10, height_50;

      for (auto it = p.begin() + 1; it != p.end(); ++it)
      {
        auto it_prev = it;
        --it_prev; //previous point
        double intensity = it->getIntensity();
        double intensity_prev = it_prev->getIntensity();
        double retention_time = it->getPos();
        double retention_time_prev = it_prev->getPos();

        // start and end intensities
        if (retention_time_prev < left && retention_time >= left)
        {
          start_intensity = intensity_prev;
        }
        else if (retention_time_prev < right && retention_time >= right)
        {
          end_intensity = intensity;
        }

        if (retention_time >= left && retention_time <= right)
        {
          //start and end retention times
          if (retention_time < peak_apex_rt_)
          {
            // start_time_at_5
            if (intensity >= 0.05*peak_height_ &&
              intensity_prev < 0.05*peak_height_ &&
              points_across_baseline_ > 1)
            {
              delta_rt = retention_time - retention_time_prev;
              delta_int = intensity - intensity_prev;
              height_5 = intensity - 0.05*peak_height_;
              start_time_at_5_ = retention_time - delta_int*delta_rt/height_5;
            }
            // start_time_at_10
            if (intensity >= 0.1*peak_height_ &&
              intensity_prev < 0.1*peak_height_ &&
              points_across_baseline_ > 1)
            {
              delta_rt = retention_time - retention_time_prev;
              delta_int = intensity - intensity_prev;
              height_10 = intensity - 0.1*peak_height_;
              start_time_at_10_ = retention_time - delta_int*delta_rt/height_10;
            }
            // start_time_at_50
            if (intensity >= 0.5*peak_height_ &&
              intensity_prev < 0.5*peak_height_ &&
              points_across_baseline_ > 1)
            {
              delta_rt = retention_time - retention_time_prev;
              delta_int = intensity - intensity_prev;
              height_50 = intensity - 0.5*peak_height_;
              start_time_at_50_ = retention_time - delta_int*delta_rt/height_50;
            }
          }
          else if (retention_time > peak_apex_rt_)
          {
            // end_time_at_5
            if (intensity <= 0.05*peak_height_ &&
              intensity_prev > 0.05*peak_height_)
            {
              delta_rt = retention_time - retention_time_prev;
              delta_int = intensity_prev - intensity;
              height_5 = 0.05*peak_height_ - intensity;
              end_time_at_5_ = retention_time - delta_int*delta_rt/height_5;
            }
            // start_time_at_10
            if (intensity <= 0.1*peak_height_ &&
              intensity_prev > 0.1*peak_height_)
            {
              delta_rt = retention_time - retention_time_prev;
              delta_int = intensity_prev - intensity;
              height_10 = 0.1*peak_height_ - intensity;
              end_time_at_10_ = retention_time - delta_int*delta_rt/height_10;
            }
            // end_time_at_50
            if (intensity <= 0.5*peak_height_ &&
            intensity_prev > 0.5*peak_height_)
            {
              delta_rt = retention_time - retention_time_prev;
              delta_int = intensity_prev - intensity;
              height_50 = 0.5*peak_height_ - intensity;
              end_time_at_50_ = retention_time - delta_int*delta_rt/height_50;
            }
          }

          // points across the peak
          points_across_baseline_++;
          if (intensity >= 0.5*peak_height_)
          {
            points_across_half_height_++;
          }
        }
      }

      // peak widths
      width_at_5_ = end_time_at_5_ - start_time_at_5_;
      width_at_10_ = end_time_at_10_ - start_time_at_10_;
      width_at_50_ = end_time_at_50_ - start_time_at_50_;
      total_width_ = right - left;
      slope_of_baseline_ = end_intensity - start_intensity;
      baseline_delta_2_height_ = slope_of_baseline_ / peak_height_;

      // other
      tailing_factor_ = width_at_5_ / std::min(peak_apex_rt_ - start_time_at_5_, end_time_at_5_ - peak_apex_rt_);
      asymmetry_factor_ = std::min(peak_apex_rt_ - start_time_at_10_, end_time_at_10_ - peak_apex_rt_) / std::max(peak_apex_rt_ - start_time_at_10_, end_time_at_10_ - peak_apex_rt_);
    }

private:
    // parameters
    String integration_type_; // intensity_sum, trapezoid, simpson
    String baseline_type_; // vertical_division, base_to_base
    String peak_model_; // none

    // outputs
    double peak_area_ = 0.0;
    double peak_height_ = -1.0;
    double peak_apex_rt_ = -1.0;
    double background_height_ = 0.0;
    double background_area_ = 0.0;
    double width_at_5_ = 0.0;
    double width_at_10_ = 0.0;
    double width_at_50_ = 0.0;
    double start_time_at_5_ = 0.0;
    double start_time_at_10_ = 0.0;
    double start_time_at_50_ = 0.0;
    double end_time_at_5_ = 0.0;
    double end_time_at_10_ = 0.0;
    double end_time_at_50_ = 0.0;
    double total_width_ = 0.0;
    /**
      The tailing factor is a measure of peak tailing.
      It is defined as the distance from the front slope of the peak to the back slope
      divided by twice the distance from the center line of the peak to the front slope,
      with all measurements made at 5% of the maximum peak height.
      tailing_factor = Tf = W0.05/2a
      where W0.05 is peak width at 5% max peak height
      a = min width to peak maximum at 5% max peak height
      b = max width to peak maximum at 5% max peak height
      0.9 < Tf < 1.2
      front Tf < 0.9
      tailing Tf > 1.2
    */
    double tailing_factor_ = 0.0;
    /**
      The asymmetry factor is a measure of peak tailing.
      It is defined as the distance from the center line of the peak to the back slope
      divided by the distance from the center line of the peak to the front slope,
      with all measurements made at 10% of the maximum peak height.
      asymmetry_factor = As = b/a
      where a is min width to peak maximum at 10% max peak height
      b is max width to peak maximum at 10% max peak height
    */
    double asymmetry_factor_ = 0.0;
    /**
      The change in baseline divided by the height is
      a way of comparing the influence of the change of baseline on the peak height.
    */
    double baseline_delta_2_height_ = 0.0;
    /**
      The slope of the baseline is a measure of slope change.
      It is approximated as the difference in baselines between the peak start and peak end.
    */
    double slope_of_baseline_ = 0.0;
    Int points_across_baseline_ = 0;
    Int points_across_half_height_ = 0;

    // helpers
    double simpson(MSChromatogram::ConstIterator it_begin, MSChromatogram::ConstIterator it_end) const;
    double simpson(MSSpectrum::ConstIterator it_begin, MSSpectrum::ConstIterator it_end) const;
  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_PEAKINTEGRATOR_H
