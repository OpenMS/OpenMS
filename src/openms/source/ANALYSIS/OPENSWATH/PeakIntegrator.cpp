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
    params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM, "The integration technique to use in integratePeak() and estimateBackground().");
    params.setValidStrings("integration_type", ListUtils::create<String>("intensity_sum,simpson,trapezoid"));
    params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE, "The baseline type to use in estimateBackground().");
    params.setValidStrings("baseline_type", ListUtils::create<String>("base_to_base,vertical_division"));
  }

  void PeakIntegrator::updateMembers_()
  {
    integration_type_ = (String)param_.getValue("integration_type");
    baseline_type_ = (String)param_.getValue("baseline_type");
  }

  template <typename PeakContainerT>
  PeakIntegrator::PeakArea PeakIntegrator::integratePeak_(const PeakContainerT& p, const double left, const double right) const
  {
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
      for (auto it = p.PosBegin(left); it != p.PosEnd(right) - 1; ++it)
      {
        peak_area += ((it + 1)->getPos() - it->getPos()) * ((it->getIntensity() + (it + 1)->getIntensity()) / 2.0);
      }
    }
    else if (integration_type_ == INTEGRATION_TYPE_SIMPSON)
    {
      if (n_points < 3)
      {
        LOG_DEBUG << std::endl << "Error in integratePeak: number of points must be >=3 for Simpson's rule" << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The number of points must be >= 3.");
      }
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
    else if (integration_type_ == INTEGRATION_TYPE_INTENSITYSUM)
    {
      std::cout << "\nWARNING: intensity_sum method is being used.\n";
      for (auto it = p.PosBegin(left); it != p.PosEnd(right); ++it)
      {
        peak_area += it->getIntensity();
      }
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
    const PeakContainerT& p, const double left, const double right,
    const double peak_apex_pos
  ) const
  {
    const double int_l = p.PosBegin(left)->getIntensity();
    const double int_r = (p.PosEnd(right) - 1)->getIntensity();
    const double delta_int = int_r - int_l;
    const double delta_pos = (p.PosEnd(right) - 1)->getPos() - p.PosBegin(left)->getPos();
    const double min_int_pos = int_r <= int_l ? (p.PosEnd(right) - 1)->getPos() : p.PosBegin(left)->getPos();
    const double delta_int_apex = std::fabs(delta_int) * std::fabs(min_int_pos - peak_apex_pos) / delta_pos;
    double background = 0.0;
    if (baseline_type_ == BASELINE_TYPE_BASETOBASE)
    {
      if (integration_type_ == INTEGRATION_TYPE_TRAPEZOID || integration_type_ == INTEGRATION_TYPE_SIMPSON)
      {
        // formula for calculating the background using the trapezoidal rule
        // background = intensity_min*delta_pos + 0.5*delta_int*delta_pos;
        background = delta_pos * (std::min(int_r, int_l) + 0.5 * std::fabs(delta_int));
      }
      else if (integration_type_ == INTEGRATION_TYPE_INTENSITYSUM)
      {
        // calculate the background using the formula
        // y = mx + b where x = rt or mz, m = slope, b = left intensity
        // sign of delta_int will determine line direction
        // background += delta_int / delta_pos * (it->getPos() - left) + int_l;
        for (auto it = p.PosBegin(left); it != p.PosEnd(right); ++it)
        {
          background += it->getPos();
        }
        UInt n_points = std::distance(p.PosBegin(left), p.PosEnd(right));
        background = (background - n_points * p.PosBegin(left)->getPos()) * delta_int / delta_pos + n_points * int_l;
      }
    }
    else if (baseline_type_ == BASELINE_TYPE_VERTICALDIVISION)
    {
      if (integration_type_ == INTEGRATION_TYPE_TRAPEZOID || integration_type_ == INTEGRATION_TYPE_SIMPSON)
      {
        background = delta_pos * std::min(int_r, int_l);
      }
      else if (integration_type_ == INTEGRATION_TYPE_INTENSITYSUM)
      {
        background = std::min(int_r, int_l) * std::distance(p.PosBegin(left), p.PosEnd(right));;
      }
    }
    else
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Please set a valid value for the parameter \"baseline_type\".");
    }
    PeakBackground pb;
    pb.area = background;
    pb.height = std::min(int_r, int_l) + delta_int_apex;
    return pb;
  }

  template <typename PeakContainerT>
  PeakIntegrator::PeakShapeMetrics PeakIntegrator::calculatePeakShapeMetrics_(
    const PeakContainerT& p, const double left, const double right,
    const double peak_height, const double peak_apex_pos
  ) const
  {
    PeakShapeMetrics psm;
    psm.points_across_baseline = 0;
    psm.points_across_half_height = 0;
    for (auto it = p.PosBegin(left); it != p.PosEnd(right); ++it)
    {
      const double intensity = it->getIntensity();
      const double intensity_prev = (it - 1)->getIntensity();
      const double position = it->getPos();
      const double position_prev = (it - 1)->getPos();
      // start and end positions (rt or mz)
      if (position < peak_apex_pos && psm.points_across_baseline > 1) // start positions
      {
        const double d_int_times_d_pos = (intensity - intensity_prev) * (position - position_prev);
        if (intensity >= 0.05 * peak_height && intensity_prev < 0.05 * peak_height)
        {
          const double height_5 = intensity - 0.05 * peak_height;
          psm.start_position_at_5 = position - d_int_times_d_pos / height_5;
        }
        if (intensity >= 0.1 * peak_height && intensity_prev < 0.1 * peak_height)
        {
          const double height_10 = intensity - 0.1 * peak_height;
          psm.start_position_at_10 = position - d_int_times_d_pos / height_10;
        }
        if (intensity >= 0.5 * peak_height && intensity_prev < 0.5 * peak_height)
        {
          const double height_50 = intensity - 0.5 * peak_height;
          psm.start_position_at_50 = position - d_int_times_d_pos / height_50;
        }
      }
      else if (position > peak_apex_pos) // end positions
      {
        const double d_int_times_d_pos = (intensity_prev - intensity) * (position - position_prev);
        if (intensity <= 0.05 * peak_height && intensity_prev > 0.05 * peak_height)
        {
          const double height_5 = 0.05 * peak_height - intensity;
          psm.end_position_at_5 = position - d_int_times_d_pos / height_5;
        }
        if (intensity <= 0.1 * peak_height && intensity_prev > 0.1 * peak_height)
        {
          const double height_10 = 0.1 * peak_height - intensity;
          psm.end_position_at_10 = position - d_int_times_d_pos / height_10;
        }
        if (intensity <= 0.5 * peak_height && intensity_prev > 0.5 * peak_height)
        {
          const double height_50 = 0.5 * peak_height - intensity;
          psm.end_position_at_50 = position - d_int_times_d_pos / height_50;
        }
      }
      // points across the peak
      ++(psm.points_across_baseline);
      if (intensity >= 0.5 * peak_height)
      {
        ++(psm.points_across_half_height);
      }
    }
    // peak widths
    psm.width_at_5 = psm.end_position_at_5 - psm.start_position_at_5;
    psm.width_at_10 = psm.end_position_at_10 - psm.start_position_at_10;
    psm.width_at_50 = psm.end_position_at_50 - psm.start_position_at_50;
    psm.total_width = (p.PosEnd(right) - 1)->getPos() - p.PosBegin(left)->getPos();
    psm.slope_of_baseline = (p.PosEnd(right) - 1)->getIntensity() - p.PosBegin(left)->getIntensity();
    psm.baseline_delta_2_height = psm.slope_of_baseline / peak_height;
    // other
    psm.tailing_factor = psm.width_at_5 / std::min(peak_apex_pos - psm.start_position_at_5, psm.end_position_at_5 - peak_apex_pos);
    psm.asymmetry_factor = std::min(peak_apex_pos - psm.start_position_at_10, psm.end_position_at_10 - peak_apex_pos) /
      std::max(peak_apex_pos - psm.start_position_at_10, psm.end_position_at_10 - peak_apex_pos);
    return psm;
  }
}
