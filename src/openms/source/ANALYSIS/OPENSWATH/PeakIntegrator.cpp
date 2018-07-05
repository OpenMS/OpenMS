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

  PeakIntegrator::PeakArea PeakIntegrator::integratePeak(const MSChromatogram& chromatogram, const double left, const double right) const
  {
    return integratePeak_(chromatogram, left, right);
  }

  PeakIntegrator::PeakArea PeakIntegrator::integratePeak(const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right) const
  {
    return integratePeak_(chromatogram, left->getRT(), right->getRT());
  }

  PeakIntegrator::PeakArea PeakIntegrator::integratePeak(const MSSpectrum& spectrum, const double left, const double right) const
  {
    return integratePeak_(spectrum, left, right);
  }

  PeakIntegrator::PeakArea PeakIntegrator::integratePeak(const MSSpectrum& spectrum, MSSpectrum::ConstIterator& left, MSSpectrum::ConstIterator& right) const
  {
    return integratePeak_(spectrum, left->getMZ(), right->getMZ());
  }

  double PeakIntegrator::simpson(MSChromatogram::ConstIterator it_begin, MSChromatogram::ConstIterator it_end) const
  {
    return simpson_(it_begin, it_end);
  }

  double PeakIntegrator::simpson(MSSpectrum::ConstIterator it_begin, MSSpectrum::ConstIterator it_end) const
  {
    return simpson_(it_begin, it_end);
  }

  PeakIntegrator::PeakBackground PeakIntegrator::estimateBackground(const MSChromatogram& chromatogram, const double left, const double right, const double peak_apex_pos) const
  {
    return estimateBackground_(chromatogram, left, right, peak_apex_pos);
  }

  PeakIntegrator::PeakBackground PeakIntegrator::estimateBackground(const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right, const double peak_apex_pos) const
  {
    return estimateBackground_(chromatogram, left->getRT(), right->getRT(), peak_apex_pos);
  }

  PeakIntegrator::PeakBackground PeakIntegrator::estimateBackground(const MSSpectrum& spectrum, const double left, const double right, const double peak_apex_pos) const
  {
    return estimateBackground_(spectrum, left, right, peak_apex_pos);
  }

  PeakIntegrator::PeakBackground PeakIntegrator::estimateBackground(const MSSpectrum& spectrum, MSSpectrum::ConstIterator& left, MSSpectrum::ConstIterator& right, const double peak_apex_pos) const
  {
    return estimateBackground_(spectrum, left->getMZ(), right->getMZ(), peak_apex_pos);
  }

  PeakIntegrator::PeakShapeMetrics PeakIntegrator::calculatePeakShapeMetrics(const MSChromatogram& chromatogram, const double left, const double right, const double peak_height, const double peak_apex_pos) const
  {
    return calculatePeakShapeMetrics_(chromatogram, left, right, peak_height, peak_apex_pos);
  }

  PeakIntegrator::PeakShapeMetrics PeakIntegrator::calculatePeakShapeMetrics(const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right, const double peak_height, const double peak_apex_pos) const
  {
    return calculatePeakShapeMetrics_(chromatogram, left->getRT(), right->getRT(), peak_height, peak_apex_pos);
  }

  PeakIntegrator::PeakShapeMetrics PeakIntegrator::calculatePeakShapeMetrics(const MSSpectrum& spectrum, const double left, const double right, const double peak_height, const double peak_apex_pos) const
  {
    return calculatePeakShapeMetrics_(spectrum, left, right, peak_height, peak_apex_pos);
  }

  PeakIntegrator::PeakShapeMetrics PeakIntegrator::calculatePeakShapeMetrics(const MSSpectrum& spectrum, MSSpectrum::ConstIterator& left, MSSpectrum::ConstIterator& right, const double peak_height, const double peak_apex_pos) const
  {
    return calculatePeakShapeMetrics_(spectrum, left->getMZ(), right->getMZ(), peak_height, peak_apex_pos);
  }

  void PeakIntegrator::getDefaultParameters(Param& params)
  {
    params.clear();

    params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM, "The integration technique to use in integratePeak() and estimateBackground() which uses either the summed intensity, integration by Simpson's rule or trapezoidal integration.");
    params.setValidStrings("integration_type", ListUtils::create<String>("intensity_sum,simpson,trapezoid"));

    params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE, "The baseline type to use in estimateBackground() based on the peak boundaries. A rectangular baseline shape is computed based either on the minimal intensity of the peak boundaries, the maximum intensity or the average intensity (base_to_base).");
    params.setValidStrings("baseline_type", ListUtils::create<String>("base_to_base,vertical_division,vertical_division_min,vertical_division_max"));

    params.setValue("fit_EMG", "false", "Fit the chromatogram/spectrum to the EMG peak model.");
    params.setValidStrings("fit_EMG", ListUtils::create<String>("false,true"));
  }

  void PeakIntegrator::updateMembers_()
  {
    integration_type_ = (String)param_.getValue("integration_type");
    baseline_type_ = (String)param_.getValue("baseline_type");
    fit_EMG_ = param_.getValue("fit_EMG").toBool();
  }

  template <typename PeakContainerT>
  PeakIntegrator::PeakArea PeakIntegrator::integratePeak_(const PeakContainerT& pc, double left, double right) const
  {
    PeakContainerT emg_pc;
    const PeakContainerT& p = EMGPreProcess_(pc, emg_pc, left, right);

    std::function<double(const double, const double)>
    compute_peak_area_trapezoid = [&p](const double left, const double right)
    {
      double peak_area { 0.0 };
      for (typename PeakContainerT::ConstIterator it = p.PosBegin(left); it != p.PosEnd(right) - 1; ++it)
      {
        peak_area += ((it + 1)->getPos() - it->getPos()) * ((it->getIntensity() + (it + 1)->getIntensity()) / 2.0);
      }
      return peak_area;
    };

    std::function<double(const double, const double)>
    compute_peak_area_intensity_sum = [&p](const double left, const double right)
    {
      // LOG_WARN << "WARNING: intensity_sum method is being used." << std::endl;
      double peak_area { 0.0 };
      for (typename PeakContainerT::ConstIterator it = p.PosBegin(left); it != p.PosEnd(right); ++it)
      {
        peak_area += it->getIntensity();
      }
      return peak_area;
    };

    double peak_area(0.0), peak_height(0.0), peak_apex_pos(0.0);
    ConvexHull2D::PointArrayType hull_points;
    UInt n_points = std::distance(p.PosBegin(left), p.PosEnd(right));
    for (auto it = p.PosBegin(left); it != p.PosEnd(right); ++it)
    {
      hull_points.push_back(DPosition<2>(it->getPos(), it->getIntensity()));
      if (peak_height < it->getIntensity())
      {
        peak_height = it->getIntensity();
        peak_apex_pos = it->getPos();
      }
    }

    if (integration_type_ == INTEGRATION_TYPE_TRAPEZOID)
    {
      if (n_points >= 2)
      {
        peak_area = compute_peak_area_trapezoid(left, right);
      }
    }
    else if (integration_type_ == INTEGRATION_TYPE_SIMPSON)
    {
      if (n_points == 2)
      {
        LOG_WARN << std::endl << "PeakIntegrator::integratePeak:"
          "number of points is 2, falling back to `trapezoid`." << std::endl;
        peak_area = compute_peak_area_trapezoid(left, right);
      }
      else if (n_points > 2)
      {
        if (n_points % 2)
        {
          peak_area = simpson(p.PosBegin(left), p.PosEnd(right));
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
              peak_area += area;
              ++valids;
            }
          }
          peak_area /= valids;
        }
      }
    }
    else if (integration_type_ == INTEGRATION_TYPE_INTENSITYSUM)
    {
      peak_area = compute_peak_area_intensity_sum(left, right);
    }
    else
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Please set a valid value for the parameter \"integration_type\".");
    }
    PeakArea pa;
    pa.area = peak_area;
    pa.height = peak_height;
    pa.apex_pos = peak_apex_pos;
    pa.hull_points = hull_points;
    return pa;
  }

  template <typename PeakContainerConstIteratorT>
  double PeakIntegrator::simpson_(PeakContainerConstIteratorT it_begin, PeakContainerConstIteratorT it_end) const
  {
    double integral = 0.0;
    for (auto it = it_begin + 1; it < it_end - 1; it = it + 2)
    {
      const double h = it->getPos() - (it - 1)->getPos();
      const double k = (it + 1)->getPos() - it->getPos();
      const double y_h = (it - 1)->getIntensity();
      const double y_0 = it->getIntensity();
      const double y_k = (it + 1)->getIntensity();
      integral += (1.0 / 6.0) * (h + k) * ((2.0 - k / h) * y_h + (pow(h + k, 2) / (h * k)) * y_0 + (2.0 - h / k) * y_k);
    }
    return integral;
  }

  template <typename PeakContainerT>
  PeakIntegrator::PeakBackground PeakIntegrator::estimateBackground_(
    const PeakContainerT& pc, double left, double right,
    const double peak_apex_pos
  ) const
  {
    PeakContainerT emg_pc;
    const PeakContainerT& p = EMGPreProcess_(pc, emg_pc, left, right);

    const double int_l = p.PosBegin(left)->getIntensity();
    const double int_r = (p.PosEnd(right) - 1)->getIntensity();
    const double delta_int = int_r - int_l;
    const double delta_pos = (p.PosEnd(right) - 1)->getPos() - p.PosBegin(left)->getPos();
    const double min_int_pos = int_r <= int_l ? (p.PosEnd(right) - 1)->getPos() : p.PosBegin(left)->getPos();
    const double delta_int_apex = std::fabs(delta_int) * std::fabs(min_int_pos - peak_apex_pos) / delta_pos;
    double area {0.0};
    double height {0.0};
    if (baseline_type_ == BASELINE_TYPE_BASETOBASE)
    {
      height = std::min(int_r, int_l) + delta_int_apex;
      if (integration_type_ == INTEGRATION_TYPE_TRAPEZOID || integration_type_ == INTEGRATION_TYPE_SIMPSON)
      {
        // formula for calculating the background using the trapezoidal rule
        // area = intensity_min*delta_pos + 0.5*delta_int*delta_pos;
        area = delta_pos * (std::min(int_r, int_l) + 0.5 * std::fabs(delta_int));
      }
      else if (integration_type_ == INTEGRATION_TYPE_INTENSITYSUM)
      {
        // calculate the background using an estimator of the form
        //    y = mx + b
        //    where x = rt or mz, m = slope, b = left intensity
        // sign of delta_int will determine line direction
        // area += delta_int / delta_pos * (it->getPos() - left) + int_l;
        double pos_sum = 0.0; // rt or mz
        for (auto it = p.PosBegin(left); it != p.PosEnd(right); ++it)
        {
          pos_sum += it->getPos();
        }
        UInt n_points = std::distance(p.PosBegin(left), p.PosEnd(right));

        // We construct the background area as the sum of a rectangular part
        // and a triangle on top. The triangle is constructed as the sum of the
        // line's y value at each sampled point: \sum_{i=0}^{n} (x_i - x_0)  * m
        const double rectangle_area = n_points * int_l;
        const double slope = delta_int / delta_pos;
        const double triangle_area = (pos_sum - n_points * p.PosBegin(left)->getPos()) * slope;
        area = triangle_area + rectangle_area;
      }
    }
    else if (baseline_type_ == BASELINE_TYPE_VERTICALDIVISION || baseline_type_ == BASELINE_TYPE_VERTICALDIVISION_MIN)
    {
      height = std::min(int_r, int_l);
      if (integration_type_ == INTEGRATION_TYPE_TRAPEZOID || integration_type_ == INTEGRATION_TYPE_SIMPSON)
      {
        area = delta_pos * std::min(int_r, int_l);
      }
      else if (integration_type_ == INTEGRATION_TYPE_INTENSITYSUM)
      {
        area = std::min(int_r, int_l) * std::distance(p.PosBegin(left), p.PosEnd(right));;
      }
    }
    else if (baseline_type_ == BASELINE_TYPE_VERTICALDIVISION_MAX)
    {
      height = std::max(int_r, int_l);
      if (integration_type_ == INTEGRATION_TYPE_TRAPEZOID || integration_type_ == INTEGRATION_TYPE_SIMPSON)
      {
        area = delta_pos * std::max(int_r, int_l);
      }
      else if (integration_type_ == INTEGRATION_TYPE_INTENSITYSUM)
      {
        area = std::max(int_r, int_l) * std::distance(p.PosBegin(left), p.PosEnd(right));
      }
    }
    else
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Please set a valid value for the parameter \"baseline_type\".");
    }
    PeakBackground pb;
    pb.area = area;
    pb.height = height;
    return pb;
  }

  template <typename PeakContainerT>
  PeakIntegrator::PeakShapeMetrics PeakIntegrator::calculatePeakShapeMetrics_(
    const PeakContainerT& pc, double left, double right,
    const double peak_height, const double peak_apex_pos
  ) const
  {
    PeakContainerT emg_pc;
    const PeakContainerT& p = EMGPreProcess_(pc, emg_pc, left, right);

    PeakShapeMetrics psm;
    psm.points_across_baseline = 0;
    psm.points_across_half_height = 0;
    for (auto it = p.PosBegin(left); it != p.PosEnd(right); ++it)
    {
      // points across the peak
      ++(psm.points_across_baseline);
      if (it->getIntensity() >= 0.5 * peak_height)
      {
        ++(psm.points_across_half_height);
      }
    }
    // positions at peak heights
    typename PeakContainerT::ConstIterator it_PosBegin_l = p.PosBegin(left);
    typename PeakContainerT::ConstIterator it_PosEnd_apex = p.PosEnd(peak_apex_pos);
    typename PeakContainerT::ConstIterator it_PosEnd_r = p.PosEnd(right);
    psm.start_position_at_5 = findPosAtPeakHeightPercent_(it_PosBegin_l, it_PosEnd_apex - 1, peak_height, 0.05, true);
    psm.start_position_at_10 = findPosAtPeakHeightPercent_(it_PosBegin_l, it_PosEnd_apex - 1, peak_height, 0.1, true);
    psm.start_position_at_50 = findPosAtPeakHeightPercent_(it_PosBegin_l, it_PosEnd_apex - 1, peak_height, 0.5, true);
    psm.end_position_at_5 = findPosAtPeakHeightPercent_(it_PosEnd_apex, it_PosEnd_r, peak_height, 0.05, false);
    psm.end_position_at_10 = findPosAtPeakHeightPercent_(it_PosEnd_apex, it_PosEnd_r, peak_height, 0.1, false);
    psm.end_position_at_50 = findPosAtPeakHeightPercent_(it_PosEnd_apex, it_PosEnd_r, peak_height, 0.5, false);
    // peak widths
    psm.width_at_5 = psm.end_position_at_5 - psm.start_position_at_5;
    psm.width_at_10 = psm.end_position_at_10 - psm.start_position_at_10;
    psm.width_at_50 = psm.end_position_at_50 - psm.start_position_at_50;
    psm.total_width = (p.PosEnd(right) - 1)->getPos() - p.PosBegin(left)->getPos();
    psm.slope_of_baseline = (p.PosEnd(right) - 1)->getIntensity() - p.PosBegin(left)->getIntensity();
    psm.baseline_delta_2_height = psm.slope_of_baseline / peak_height;
    // Source of tailing_factor and asymmetry_factor formulas:
    // USP 40 - NF 35 The United States Pharmacopeia and National Formulary - Supplementary
    psm.tailing_factor = psm.width_at_5 / (2*(peak_apex_pos - psm.start_position_at_5));
    psm.asymmetry_factor = (psm.end_position_at_10 - peak_apex_pos) / (peak_apex_pos - psm.start_position_at_10);
    return psm;
  }

  template <typename PeakContainerConstIteratorT>
  double PeakIntegrator::findPosAtPeakHeightPercent_(
    PeakContainerConstIteratorT it_begin,
    PeakContainerConstIteratorT it_end,
    const double peak_height,
    const double percent,
    const bool is_left_half
  ) const
  {
    const double percent_intensity = peak_height * percent;
    PeakContainerConstIteratorT closest = is_left_half ? it_begin : it_end - 1;
    if (is_left_half)
    {
      for (
        PeakContainerConstIteratorT it = it_begin;
        it != it_end && it->getIntensity() <= percent_intensity;
        closest = it++
      ) {}
    }
    else
    {
      for (
        PeakContainerConstIteratorT it = it_end - 1;
        it != it_begin && it->getIntensity() <= percent_intensity;
        closest = it--
      ) {}
    }
    return closest->getPos();
  }

  template <typename PeakContainerT>
  const PeakContainerT& PeakIntegrator::EMGPreProcess_(
    const PeakContainerT& pc,
    PeakContainerT& emg_pc,
    double& left,
    double& right
  ) const
  {
    if (fit_EMG_)
    {
      emg_.fitEMGPeakModel(pc, emg_pc, left, right);
      left = emg_pc.front().getPos();
      right = emg_pc.back().getPos();
      return emg_pc;
    }
    return pc;
  }
}
