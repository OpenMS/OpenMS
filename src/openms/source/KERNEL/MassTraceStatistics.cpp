// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MassTraceStatistics.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace OpenMS::MassTraceStatistics
{
  namespace
  {
    /// Median of @p values (which get sorted). @pre values is not empty
    double medianOf(std::vector<double>& values)
    {
      std::sort(values.begin(), values.end());
      const Size mid = values.size() / 2;
      return (values.size() % 2 == 0) ? (values[mid - 1] + values[mid]) / 2 : values[mid];
    }

    /// Index of the first maximum of intensity(0), ..., intensity(size - 1). @pre size > 0
    template <typename IntensityAt>
    Size indexOfMax(Size size, const IntensityAt& intensity)
    {
      double max_int = intensity(0);
      Size max_idx(0);
      for (Size i = 1; i < size; ++i)
      {
        const double act_int = intensity(i);
        if (act_int > max_int)
        {
          max_int = act_int;
          max_idx = i;
        }
      }
      return max_idx;
    }

    /// Area under intensity(i) over the RTs of @p peaks (trapezoidal rule)
    template <typename IntensityAt>
    double trapezoidArea(PeakRange peaks, const IntensityAt& intensity)
    {
      double area(0.0);
      for (Size i = 1; i < peaks.size(); ++i)
      {
        area += (intensity(i - 1) + intensity(i)) / 2 * (peaks[i].getRT() - peaks[i - 1].getRT());
      }
      return area;
    }

    /// calculate x coordinate of start/end indexes at half_max
    /// calculation is based on (yB - yA) / (xB - xA) = (y_eval - yA) / (xC - xA)
    /// solve for xC: xC = xA + ((y_eval - yA) * (xB - xA) / (yB - yA))
    /// (xA, yA) and (xB, yB) must be corresponding points of the same line, i.e.
    /// each x has to be passed together with its own y. y must ascend from yA to
    /// yB (see the precondition), x may ascend or descend; xC is returned between
    /// xA and xB either way.
    double linearInterpolationAtY(double xA, double xB, double yA, double yB, double y_eval)
    {
      OPENMS_PRECONDITION(yA <= y_eval && y_eval <= yB, "y_eval is not between yA and yB")
      // no solution -> return an estimate
      if (std::fabs(xA - xB) == 0 || std::fabs(yA - yB) == 0)  { return xA; }

      double xC = (xA + ((y_eval - yA) * (xB - xA) / (yB - yA)));
      // the interpolation parameter is in [0, 1] because of the precondition, so xC lies
      // between xA and xB -- in either x order, since x need not ascend with y
      OPENMS_POSTCONDITION(std::min(xA, xB) <= xC && xC <= std::max(xA, xB), "xC is not between xA and xB");

      return xC;
    }

    /// FWHM around the peak @p max_idx, based on intensity(0), ..., intensity(size - 1). @pre size > 0
    template <typename IntensityAt>
    FWHM estimateFWHMAround(PeakRange peaks, Size max_idx, Size size, const IntensityAt& intensity)
    {
      // The flank search below walks peak indices outward from the apex, so the trace
      // has to be sorted along its x axis. This used to be caught only indirectly, by
      // the postcondition of linearInterpolationAtY firing when a shuffled axis
      // happened to hand it a descending bracketing pair (see #10051, where a
      // non-stable sort on equal RT permutes the ion-mobility axis PeakPickerIM builds
      // its traces from). That postcondition is not a sortedness check -- it is correct
      // for either x order -- so the invariant is stated here instead, where it holds
      // for every shuffle rather than only the ones that invert a bracket.
      OPENMS_PRECONDITION(std::is_sorted(peaks.begin(), peaks.end(),
                                         [](const Peak2D& a, const Peak2D& b) { return a.getRT() < b.getRT(); }),
                          "mass trace is not sorted along its x (RT) axis")

      FWHM result;

      // no points left of apex in mass trace OR no points right of apex in mass trace
      if (max_idx == 0 || max_idx == size - 1)
      {
        return result;
      }

      const double half_max_int(intensity(max_idx) / 2.0);

      Size left_border(max_idx), right_border(max_idx);

      while (left_border > 0 && intensity(left_border) >= half_max_int)
      {
        --left_border;
      }

      while (right_border + 1 < size && intensity(right_border) >= half_max_int)
      {
        ++right_border;
      }

      result.start_idx = left_border;
      result.end_idx = right_border;

      // calculate RT of left side
      // the pair below half max is (left_rt_bottom, intensity(left_border)),
      // the pair above half max is (left_rt_top, intensity(fwhm_left_top_idx))
      Size fwhm_left_top_idx = left_border + 1;
      double left_rt_bottom = peaks[left_border].getRT();
      double left_rt_top = peaks[fwhm_left_top_idx].getRT();
      double fwhm_begin_rt;
      // no points below half max -> don't extrapolate, use leftmost RT
      if (intensity(left_border) > half_max_int)
      {
        fwhm_begin_rt = peaks[left_border].getRT();
      }
      else
      {
        fwhm_begin_rt = linearInterpolationAtY(left_rt_bottom, left_rt_top,
                                               intensity(left_border),
                                               intensity(fwhm_left_top_idx),
                                               half_max_int);
      }

      // calculate RT of right side
      // as on the left side, the point below half max comes first: it is
      // (right_rt_bottom, intensity(right_border)), the one above half max is
      // (right_rt_top, intensity(fwhm_right_top_idx)). Each x must travel with
      // its own y, so the x arguments are in descending order here for a trace
      // with ascending RT.
      Size fwhm_right_top_idx = right_border - 1;
      double right_rt_bottom = peaks[right_border].getRT();
      double right_rt_top = peaks[fwhm_right_top_idx].getRT();

      double fwhm_end_rt;
      // no points below half max -> don't extrapolate, use rightmost RT
      if (intensity(right_border) > half_max_int)
      {
        fwhm_end_rt = peaks[right_border].getRT();
      }
      else
      {
        fwhm_end_rt = linearInterpolationAtY(right_rt_bottom, right_rt_top,
                                             intensity(right_border),
                                             intensity(fwhm_right_top_idx),
                                             half_max_int);
      }

      result.fwhm = std::fabs(fwhm_end_rt - fwhm_begin_rt);

      OPENMS_POSTCONDITION(result.fwhm <= std::fabs(peaks[right_border].getRT() - peaks[left_border].getRT()), "Interpolated fwhm is not smaller than extrapolated fwhm!")

      return result;
    }
  } // anonymous namespace

  double meanMZ(PeakRange peaks)
  {
    if (peaks.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace is empty... centroid MZ undefined!", StringUtils::toStr(peaks.size()));
    }

    if (peaks.size() == 1)
    {
      return peaks.front().getMZ();
    }

    double sum_mz(0.0);
    for (const Peak2D& peak : peaks)
    {
      sum_mz += peak.getMZ();
    }

    return sum_mz / peaks.size();
  }

  double medianMZ(PeakRange peaks)
  {
    if (peaks.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace is empty... centroid MZ undefined!", StringUtils::toStr(peaks.size()));
    }

    if (peaks.size() == 1)
    {
      return peaks.front().getMZ();
    }

    std::vector<double> mzs;
    mzs.reserve(peaks.size());
    for (const Peak2D& peak : peaks)
    {
      mzs.push_back(peak.getMZ());
    }

    return medianOf(mzs);
  }

  double weightedMeanMZ(PeakRange peaks)
  {
    if (peaks.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace is empty... centroid MZ undefined!", StringUtils::toStr(peaks.size()));
    }

    if (peaks.size() == 1)
    {
      return peaks.front().getMZ();
    }

    double weighted_sum(0.0);
    double total_weight(0.0);

    for (const Peak2D& peak : peaks)
    {
      double w_i = peak.getIntensity();
      total_weight += w_i;
      weighted_sum += w_i * peak.getMZ();
    }

    if (total_weight < std::numeric_limits<double>::epsilon())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "All weights were equal to zero! Empty trace? Aborting...", StringUtils::toStr(total_weight));
    }

    return weighted_sum / total_weight;
  }

  double weightedMZsd(PeakRange peaks, double centroid_mz)
  {
    if (peaks.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace is empty... std of MZ undefined!", StringUtils::toStr(peaks.size()));
    }

    double weighted_sum(0.0);
    double total_weight(0.0);

    for (const Peak2D& peak : peaks)
    {
      double w_i = peak.getIntensity();
      total_weight += w_i;
      weighted_sum += w_i * std::exp(2 * std::log(std::abs(peak.getMZ() - centroid_mz)));
    }

    if (total_weight < std::numeric_limits<double>::epsilon())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "All weights were equal to zero! Empty trace? Aborting...", StringUtils::toStr(total_weight));
    }

    return std::sqrt(weighted_sum) / std::sqrt(total_weight);
  }

  double weightedMeanRT(PeakRange peaks)
  {
    if (peaks.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace is empty... centroid RT undefined!", StringUtils::toStr(peaks.size()));
    }

    // support single peak traces
    if (peaks.size() == 1)
    {
      return peaks.front().getRT();
    }

    // Intensity-weighted mean RT, integrated with the trapezoidal rule: each peak is
    // weighted by its intensity times the RT span it stands for, which is half the
    // distance to each of its neighbours (the outermost peaks only have the half towards
    // their single neighbour).
    //
    // The previous formulation weighted every peak by the distance to its *predecessor*
    // and started the sum at the second peak, so peak 0 contributed to neither the
    // numerator nor the denominator. That biased the centroid towards later RT and, for a
    // two-point trace, degenerated to the RT of the second peak whatever the intensities
    // (see issue #2777).
    double wmean_rt(0.0);
    double trace_area(0.0);

    for (Size i = 0; i < peaks.size(); ++i)
    {
      const double rt_lower = (i > 0) ? peaks[i - 1].getRT() : peaks[i].getRT();
      const double rt_upper = (i + 1 < peaks.size()) ? peaks[i + 1].getRT() : peaks[i].getRT();
      const double weight = peaks[i].getIntensity() * (rt_upper - rt_lower) / 2.0;
      wmean_rt += weight * peaks[i].getRT();
      trace_area += weight;
    }

    if (trace_area < std::numeric_limits<double>::epsilon())
    { // no intensity at all, or all peaks share one RT -> weights carry no information
      double rt_sum(0.0);
      for (const Peak2D& peak : peaks)
      {
        rt_sum += peak.getRT();
      }
      return rt_sum / peaks.size();
    }

    return wmean_rt / trace_area;
  }

  double smoothedWeightedMeanRT(PeakRange peaks, IntensityRange smoothed)
  {
    if (smoothed.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace was not smoothed before! Aborting...", StringUtils::toStr(smoothed.size()));
    }
    OPENMS_PRECONDITION(smoothed.size() == peaks.size(), "Number of smoothed intensities deviates from the number of peaks")

    // support single peak traces
    if (peaks.size() == 1)
    {
      return peaks.front().getRT();
    }

    double trace_area(0.0), wmean_rt(0.0);

    for (Size i = 0; i < smoothed.size(); ++i)
    {
      if (smoothed[i] > 0.0)
      {
        wmean_rt += smoothed[i] * peaks[i].getRT();
        trace_area += smoothed[i];
      }
    }

    if (trace_area < std::numeric_limits<double>::epsilon())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peak area equals to zero... impossible to compute weights!", StringUtils::toStr(peaks.size()));
    }

    return wmean_rt / trace_area;
  }

  double smoothedMaxRT(PeakRange peaks, IntensityRange smoothed)
  {
    if (smoothed.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace was not smoothed before! Aborting...", StringUtils::toStr(smoothed.size()));
    }
    OPENMS_PRECONDITION(smoothed.size() == peaks.size(), "Number of smoothed intensities deviates from the number of peaks")

    // Handle single-peak traces: the centroid RT is the RT of the only peak
    if (peaks.size() == 1)
    {
      return peaks.front().getRT();
    }

    double tmp_max(-1.0);
    Size max_idx(0);

    for (Size i = 0; i < smoothed.size(); ++i)
    {
      if (smoothed[i] > tmp_max)
      {
        tmp_max = smoothed[i];
        max_idx = i;
      }
    }

    if (tmp_max <= 0.0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Negative max intensity encountered!", StringUtils::toStr(tmp_max));
    }

    return peaks[max_idx].getRT();
  }

  double medianRT(PeakRange peaks)
  {
    if (peaks.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace is empty... centroid RT undefined!", StringUtils::toStr(peaks.size()));
    }

    if (peaks.size() == 1)
    {
      return peaks.front().getRT();
    }

    std::vector<double> rts;
    rts.reserve(peaks.size());
    for (const Peak2D& peak : peaks)
    {
      rts.push_back(peak.getRT());
    }

    return medianOf(rts);
  }

  Size apexIndex(PeakRange peaks)
  {
    if (peaks.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace appears to be empty! Aborting...", StringUtils::toStr(peaks.size()));
    }

    return indexOfMax(peaks.size(), [&peaks](Size i) -> double { return peaks[i].getIntensity(); });
  }

  Size apexIndex(PeakRange peaks, IntensityRange smoothed)
  {
    if (smoothed.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace was not smoothed before! Aborting...", StringUtils::toStr(smoothed.size()));
    }

    if (peaks.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace appears to be empty! Aborting...", StringUtils::toStr(peaks.size()));
    }
    OPENMS_PRECONDITION(smoothed.size() == peaks.size(), "Number of smoothed intensities deviates from the number of peaks")

    return indexOfMax(peaks.size(), [&smoothed](Size i) { return smoothed[i]; });
  }

  FWHM estimateFWHM(PeakRange peaks)
  {
    const Size max_idx = apexIndex(peaks);
    return estimateFWHMAround(peaks, max_idx, peaks.size(), [&peaks](Size i) -> double { return peaks[i].getIntensity(); });
  }

  FWHM estimateFWHM(PeakRange peaks, IntensityRange smoothed)
  {
    const Size max_idx = apexIndex(peaks, smoothed);
    return estimateFWHMAround(peaks, max_idx, smoothed.size(), [&smoothed](Size i) { return smoothed[i]; });
  }

  double intensitySum(PeakRange peaks)
  {
    double sum(0.0);
    for (const Peak2D& peak : peaks)
    {
      sum += peak.getIntensity();
    }
    return sum;
  }

  double maxIntensity(PeakRange peaks)
  {
    double max_int(0.0);
    for (const Peak2D& peak : peaks)
    {
      if (peak.getIntensity() > max_int)
      {
        max_int = peak.getIntensity();
      }
    }
    return max_int;
  }

  double maxIntensity(IntensityRange intensities)
  {
    double max_int(0.0);
    for (const double intensity : intensities)
    {
      if (intensity > max_int)
      {
        max_int = intensity;
      }
    }
    return max_int;
  }

  double medianIntensity(PeakRange peaks)
  {
    OPENMS_PRECONDITION(!peaks.empty(), "Median intensity of an empty mass trace is undefined")

    std::vector<double> ints;
    ints.reserve(peaks.size());
    for (const Peak2D& peak : peaks)
    {
      ints.push_back(peak.getIntensity());
    }

    return medianOf(ints);
  }

  double peakArea(PeakRange peaks)
  {
    return trapezoidArea(peaks, [&peaks](Size i) -> double { return peaks[i].getIntensity(); });
  }

  double peakArea(PeakRange peaks, IntensityRange smoothed)
  {
    OPENMS_PRECONDITION(smoothed.size() == peaks.size(), "Number of smoothed intensities deviates from the number of peaks")
    return trapezoidArea(peaks, [&smoothed](Size i) { return smoothed[i]; });
  }
} // namespace OpenMS::MassTraceStatistics
