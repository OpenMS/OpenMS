// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/MATH/MISC/EmgGradientDescent.h>

namespace OpenMS
{
  /**
    @brief Compute the area, background and shape metrics of a peak.

    The area computation is performed in integratePeak() and it supports
    integration by simple sum of the intensity, integration by Simpson's rule
    implementations for an odd number of unequally spaced points or integration
    by the trapezoid rule.

    The background computation is performed in estimateBackground() and it
    supports three different approaches to baseline correction, namely
    computing a rectangular shape under the peak based on the minimum value of
    the peak borders (vertical_division_min), a rectangular shape based on the
    maximum value of the beak borders (vertical_division_max) or a trapezoidal
    shape based on a straight line between the peak borders (base_to_base).

    Peak shape metrics are computed in calculatePeakShapeMetrics() and multiple
    metrics are supported.

    The containers supported by the methods are MSChromatogram and MSSpectrum.
    
    @warning integratePeak() using Simpson's rule can result in negative areas despite
    strictly positive intensities in the input dataset. An example is given in
    the class test (see area = -665788.77). 
  */
  class OPENMS_DLLAPI PeakIntegrator :
    public DefaultParamHandler
  {
public:
    /// Constructor
    PeakIntegrator();
    /// Destructor
    ~PeakIntegrator() override;

    /** @name integratePeak() output
      The integratePeak() method uses this struct to save its results.
    */
    ///@{
    struct PeakArea
    {
      /**
        The peak's computed area
      */
      double area = 0.0;
      /**
        The peak's highest intensity
      */
      double height = 0.0;
      /**
        The position of the point with highest intensity
      */
      double apex_pos = 0.0;
      /**
        The peak's hull points
      */
      ConvexHull2D::PointArrayType hull_points;
    };
    ///@}

    /** @name estimateBackground() output
      The estimateBackground() method uses this struct to save its results.
    */
    ///@{
    struct PeakBackground
    {
      /**
        The background area estimation
      */
      double area = 0.0;
      /**
        The background height
      */
      double height = 0.0;
    };
    ///@}

    /** @name calculatePeakShapeMetrics() output
      
        The calculatePeakShapeMetrics() method uses this struct to save its results.
    */
    ///@{
    struct PeakShapeMetrics
    {
      /**
        The width of the peak at 5% the peak's height.
      */
      double width_at_5 = 0.0;
      /**
        The width of the peak at 10% the peak's height.
      */
      double width_at_10 = 0.0;
      /**
        The width of the peak at 50% the peak's height.
      */
      double width_at_50 = 0.0;
      /**
        The start position at which the intensity is 5% the peak's height.
      */
      double start_position_at_5 = 0.0;
      /**
        The start position at which the intensity is 10% the peak's height.
      */
      double start_position_at_10 = 0.0;
      /**
        The start position at which the intensity is 50% the peak's height.
      */
      double start_position_at_50 = 0.0;
      /**
        The end position at which the intensity is 5% the peak's height.
      */
      double end_position_at_5 = 0.0;
      /**
        The end position at which the intensity is 10% the peak's height.
      */
      double end_position_at_10 = 0.0;
      /**
        The end position at which the intensity is 50% the peak's height.
      */
      double end_position_at_50 = 0.0;
      /**
        The peak's total width.
      */
      double total_width = 0.0;
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
      double tailing_factor = 0.0;
      /**
        The asymmetry factor is a measure of peak tailing.
        It is defined as the distance from the center line of the peak to the back slope
        divided by the distance from the center line of the peak to the front slope,
        with all measurements made at 10% of the maximum peak height.
        asymmetry_factor = As = b/a
        where a is min width to peak maximum at 10% max peak height
        b is max width to peak maximum at 10% max peak height
      */
      double asymmetry_factor = 0.0;
      /**
        The slope of the baseline is a measure of slope change.
        It is approximated as the difference in baselines between the peak start and peak end.
      */
      double slope_of_baseline = 0.0;
      /**
        The change in baseline divided by the height is
        a way of comparing the influence of the change of baseline on the peak height.
      */
      double baseline_delta_2_height = 0.0;
      /**
        The number of points across the baseline.
      */
      Int points_across_baseline = 0;
      /**
        The number of points across half the peak's height.
      */
      Int points_across_half_height = 0;
    };
    ///@}

    /** @name Constant expressions for parameters
      
        Constants expressions used throughout the code and tests to set
        the integration and baseline types.
    */
    ///@{
    /// Integration type: intensity sum
    static constexpr const char* INTEGRATION_TYPE_INTENSITYSUM = "intensity_sum";
    /// Integration type: trapezoid
    static constexpr const char* INTEGRATION_TYPE_TRAPEZOID = "trapezoid";
    /// Integration type: simpson
    static constexpr const char* INTEGRATION_TYPE_SIMPSON = "simpson";
    /// Baseline type: base to base
    static constexpr const char* BASELINE_TYPE_BASETOBASE = "base_to_base";
    /// Baseline type: vertical division (min of end points; only for backwards compatibility)
    static constexpr const char* BASELINE_TYPE_VERTICALDIVISION = "vertical_division";
    /// Baseline type: vertical division (min of end points)
    static constexpr const char* BASELINE_TYPE_VERTICALDIVISION_MIN = "vertical_division_min";
    /// Baseline type: vertical division (max of end points)
    static constexpr const char* BASELINE_TYPE_VERTICALDIVISION_MAX = "vertical_division_max";
    ///@}

    /**
      @brief Compute the area of a peak contained in a MSChromatogram.

      The value of integration_type_ decides which integration technique to use:
      - "trapezoid" for the trapezoidal rule
      - "simpson" for the Simpson's rule (for unequally spaced points, Shklov, 1960)
      - "intensity_sum" for the simple sum of the intensities

      @note Make sure the chromatogram is sorted with respect to retention time.

      @throw Exception::InvalidParameter for class parameter `integration_type`.

      @param[in] chromatogram The chromatogram which contains the peak
      @param[in] left The left retention time boundary
      @param[in] right The right retention time boundary

      @return A struct containing the informations about the peak's area, height and position
    */
    PeakArea integratePeak(
      const MSChromatogram& chromatogram, const double left, const double right
    ) const;

    /**
      @brief Compute the area of a peak contained in a MSChromatogram.

      The value of integration_type_ decides which integration technique to use:
      - "trapezoid" for the trapezoidal rule
      - "simpson" for the Simpson's rule (for unequally spaced points, Shklov, 1960)
      - "intensity_sum" for the simple sum of the intensities

      @note Make sure the chromatogram is sorted with respect to retention time.

      @throw Exception::InvalidParameter for class parameter `integration_type`.

      @param[in] chromatogram The chromatogram which contains the peak
      @param[in] left The iterator to the first point
      @param[in] right The iterator to the last point

      @return A struct containing the informations about the peak's area, height and position
    */
    PeakArea integratePeak(
      const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right
    ) const;

    /**
      @brief Compute the area of a peak contained in a MSSpectrum.

      The value of integration_type_ decides which integration technique to use:
      - "trapezoid" for the trapezoidal rule
      - "simpson" for the Simpson's rule (for unequally spaced points, Shklov, 1960)
      - "intensity_sum" for the simple sum of the intensities

      @note Make sure the spectrum is sorted with respect to mass-to-charge ratio.

      @throw Exception::InvalidParameter for class parameter `integration_type`.

      @param[in] spectrum The spectrum which contains the peak
      @param[in] left The left mass-to-charge ratio boundary
      @param[in] right The right mass-to-charge ratio boundary

      @return A struct containing the informations about the peak's area, height and position
    */
    PeakArea integratePeak(
      const MSSpectrum& spectrum, const double left, const double right
    ) const;

    /**
      @brief Compute the area of a peak contained in a MSSpectrum.

      The value of integration_type_ decides which integration technique to use:
      - "trapezoid" for the trapezoidal rule
      - "simpson" for the Simpson's rule (for unequally spaced points, Shklov, 1960)
      - "intensity_sum" for the simple sum of the intensities

      @note Make sure the spectrum is sorted with respect to mass-to-charge ratio.

      @throw Exception::InvalidParameter for class parameter `integration_type`.

      @param[in] spectrum The spectrum which contains the peak
      @param[in] left The iterator to the first point
      @param[in] right The iterator to the last point

      @return A struct containing the informations about the peak's area, height and position
    */
    PeakArea integratePeak(
      const MSSpectrum& spectrum, MSSpectrum::ConstIterator& left, MSSpectrum::ConstIterator& right
    ) const;

    /**
      @brief Estimate the background of a peak contained in a MSChromatogram.

      The user can choose to compute one of two background types: "vertical_sum" and "base_to_base".
      For the former case, the area is computed as a rectangle with delta RT being the base and
      the minimum intensity on boundaries as the height.
      For the latter case, the area is computed as a rectangle trapezoid. Similar to the "vertical_sum"
      solution, this technique also takes into account the area between the intensities on boundaries.

      For both cases, the parameter integration_type_ decides which formula to use to compute the area.
      The user should make sure to use the same integration_type between calls of estimateBackground() and
      integratePeak().

      @note Make sure the chromatogram is sorted with respect to retention time.

      @throw Exception::InvalidParameter for class parameter `baseline_type`.

      @param[in] chromatogram The chromatogram which contains the peak
      @param[in] left The left retention time boundary
      @param[in] right The right retention time boundary
      @param[in] peak_apex_pos The position of the point with highest intensity

      @return A struct containing the informations about the peak's background area and height
    */
    PeakBackground estimateBackground(
      const MSChromatogram& chromatogram, const double left, const double right,
      const double peak_apex_pos
    ) const;

    /**
      @brief Estimate the background of a peak contained in a MSChromatogram.

      The user can choose to compute one of two background types: "vertical_sum" and "base_to_base".
      For the former case, the area is computed as a rectangle with delta RT being the base and
      the minimum intensity on boundaries as the height.
      For the latter case, the area is computed as a rectangle trapezoid. Similar to the "vertical_sum"
      solution, this technique also takes into account the area between the intensities on boundaries.

      For both cases, the parameter integration_type_ decides which formula to use to compute the area.
      The user should make sure to use the same integration_type between calls of estimateBackground() and
      integratePeak().

      @note Make sure the chromatogram is sorted with respect to retention time.

      @throw Exception::InvalidParameter for class parameter `baseline_type`.

      @param[in] chromatogram The chromatogram which contains the peak
      @param[in] left The iterator to the first point
      @param[in] right The iterator to the last point
      @param[in] peak_apex_pos The position of the point with highest intensity

      @return A struct containing the informations about the peak's background area and height
    */
    PeakBackground estimateBackground(
      const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right,
      const double peak_apex_pos
    ) const;

    /**
      @brief Estimate the background of a peak contained in a MSSpectrum.

      The user can choose to compute one of two background types: "vertical_sum" and "base_to_base".
      For the former case, the area is computed as a rectangle with delta MZ being the base and
      the minimum intensity on boundaries as the height.
      For the latter case, the area is computed as a rectangle trapezoid. Similar to the "vertical_sum"
      solution, this technique also takes into account the area between the intensities on boundaries.

      For both cases, the parameter integration_type_ decides which formula to use to compute the area.
      The user should make sure to use the same integration_type between calls of estimateBackground() and
      integratePeak().

      @note Make sure the spectrum is sorted with respect to mass-to-charge ratio.

      @throw Exception::InvalidParameter for class parameter `baseline_type`.

      @param[in] spectrum The spectrum which contains the peak
      @param[in] left The left mass-to-charge ratio boundary
      @param[in] right The right mass-to-charge ratio boundary
      @param[in] peak_apex_pos The position of the point with highest intensity

      @return A struct containing the informations about the peak's background area and height
    */
    PeakBackground estimateBackground(
      const MSSpectrum& spectrum, const double left, const double right,
      const double peak_apex_pos
    ) const;

    /**
      @brief Estimate the background of a peak contained in a MSSpectrum.

      The user can choose to compute one of two background types: "vertical_sum" and "base_to_base".
      For the former case, the area is computed as a rectangle with delta MZ being the base and
      the minimum intensity on boundaries as the height.
      For the latter case, the area is computed as a rectangle trapezoid. Similar to the "vertical_sum"
      solution, this technique also takes into account the area between the intensities on boundaries.

      For both cases, the parameter integration_type_ decides which formula to use to compute the area.
      The user should make sure to use the same integration_type between calls of estimateBackground() and
      integratePeak().

      @note Make sure the spectrum is sorted with respect to mass-to-charge ratio.

      @throw Exception::InvalidParameter for class parameter `baseline_type`.

      @param[in] spectrum The spectrum which contains the peak
      @param[in] left The iterator to the first point
      @param[in] right The iterator to the last point
      @param[in] peak_apex_pos The position of the point with highest intensity

      @return A struct containing the informations about the peak's background area and height
    */
    PeakBackground estimateBackground(
      const MSSpectrum& spectrum, MSSpectrum::ConstIterator& left, MSSpectrum::ConstIterator& right,
      const double peak_apex_pos
    ) const;

    /**
      @brief Calculate peak's shape metrics.

      The calculated characteristics are the start and end times at 0.05, 0.10 and
      0.5 the peak's height. Also the widths at those positions are calculated.
      Other values: the peak's total width, its tailing factor, asymmetry factor,
      baseline delta to height and the slope of the baseline.
      The number of points across the baseline and also at half height are saved.

      @note Make sure the chromatogram is sorted with respect to retention time.

      @param[in] chromatogram The chromatogram which contains the peak
      @param[in] left The left retention time boundary
      @param[in] right The right retention time boundary
      @param[in] peak_height The peak's highest intensity
      @param[in] peak_apex_pos The position of the point with highest intensity

      @return A struct containing the calculated peak shape metrics
    */
    PeakShapeMetrics calculatePeakShapeMetrics(
      const MSChromatogram& chromatogram, const double left, const double right,
      const double peak_height, const double peak_apex_pos
    ) const;

    /**
      @brief Calculate peak's shape metrics.

      The calculated characteristics are the start and end times at 0.05, 0.10 and
      0.5 the peak's height. Also the widths at those positions are calculated.
      Other values: the peak's total width, its tailing factor, asymmetry factor,
      baseline delta to height and the slope of the baseline.
      The number of points across the baseline and also at half height are saved.

      @note Make sure the chromatogram is sorted with respect to retention time.

      @param[in] chromatogram The chromatogram which contains the peak
      @param[in] left The iterator to the first point
      @param[in] right The iterator to the last point
      @param[in] peak_height The peak's highest intensity
      @param[in] peak_apex_pos The position of the point with highest intensity

      @return A struct containing the calculated peak shape metrics
    */
    PeakShapeMetrics calculatePeakShapeMetrics(
      const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right,
      const double peak_height, const double peak_apex_pos
    ) const;

    /**
      @brief Calculate peak's shape metrics.

      The calculated characteristics are the start and end positions at 0.05, 0.10 and
      0.5 the peak's height. Also the widths at those positions are calculated.
      Other values: the peak's total width, its tailing factor, asymmetry factor,
      baseline delta to height and the slope of the baseline.
      The number of points across the baseline and also at half height are saved.

      @note Make sure the spectrum is sorted with respect to mass-to-charge ratio.

      @param[in] spectrum The spectrum which contains the peak
      @param[in] left The left mass-to-charge ratio boundary
      @param[in] right The right mass-to-charge ratio boundary
      @param[in] peak_height The peak's highest intensity
      @param[in] peak_apex_pos The position of the point with highest intensity

      @return A struct containing the calculated peak shape metrics
    */
    PeakShapeMetrics calculatePeakShapeMetrics(
      const MSSpectrum& spectrum, const double left, const double right,
      const double peak_height, const double peak_apex_pos
    ) const;

    /**
      @brief Calculate peak's shape metrics.

      The calculated characteristics are the start and end positions at 0.05, 0.10 and
      0.5 the peak's height. Also the widths at those positions are calculated.
      Other values: the peak's total width, its tailing factor, asymmetry factor,
      baseline delta to height and the slope of the baseline.
      The number of points across the baseline and also at half height are saved.

      @note Make sure the spectrum is sorted with respect to mass-to-charge ratio.

      @param[in] spectrum The spectrum which contains the peak
      @param[in] left The iterator to the first point
      @param[in] right The iterator to the last point
      @param[in] peak_height The peak's highest intensity
      @param[in] peak_apex_pos The position of the point with highest intensity

      @return A struct containing the calculated peak shape metrics
    */
    PeakShapeMetrics calculatePeakShapeMetrics(
      const MSSpectrum& spectrum, MSSpectrum::ConstIterator& left, MSSpectrum::ConstIterator& right,
      const double peak_height, const double peak_apex_pos
    ) const;

    void getDefaultParameters(Param& params);

protected:
    void updateMembers_() override;

    template <typename PeakContainerT>
    PeakArea integratePeak_(const PeakContainerT& pc, double left, double right) const
    {
      OPENMS_PRECONDITION(left <= right, "Left peak boundary must be smaller than right boundary!") // otherwise the code below will segfault (due to PosBegin/PosEnd)
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
        // OPENMS_LOG_WARN << "WARNING: intensity_sum method is being used." << std::endl;
        double peak_area { 0.0 };
        for (typename PeakContainerT::ConstIterator it = p.PosBegin(left); it != p.PosEnd(right); ++it)
        {
          peak_area += it->getIntensity();
        }
        return peak_area;
      };

      PeakArea pa;
      pa.apex_pos = (left + right) / 2; // initial estimate, to avoid apex being outside of [left,right]
      UInt n_points = std::distance(p.PosBegin(left), p.PosEnd(right));
      for (auto it = p.PosBegin(left); it != p.PosEnd(right); ++it) //OMS_CODING_TEST_EXCLUDE
      {
        pa.hull_points.push_back(DPosition<2>(it->getPos(), it->getIntensity()));
        if (pa.height < it->getIntensity())
        {
          pa.height = it->getIntensity();
          pa.apex_pos = it->getPos();
        }
      }

      if (integration_type_ == INTEGRATION_TYPE_TRAPEZOID)
      {
        if (n_points >= 2)
        {
          pa.area = compute_peak_area_trapezoid(left, right);
        }
      }
      else if (integration_type_ == INTEGRATION_TYPE_SIMPSON)
      {
        if (n_points == 2)
        {
          OPENMS_LOG_WARN << std::endl << "PeakIntegrator::integratePeak:"
            "number of points is 2, falling back to `trapezoid`." << std::endl;
          pa.area = compute_peak_area_trapezoid(left, right);
        }
        else if (n_points > 2)
        {
          if (n_points % 2)
          {
            pa.area = simpson_(p.PosBegin(left), p.PosEnd(right));
          }
          else
          {
            double areas[4] = {-1.0, -1.0, -1.0, -1.0};
            areas[0] = simpson_(p.PosBegin(left), p.PosEnd(right) - 1);   // without last point
            areas[1] = simpson_(p.PosBegin(left) + 1, p.PosEnd(right));   // without first point
            if (p.begin() <= p.PosBegin(left) - 1)
            {
              areas[2] = simpson_(p.PosBegin(left) - 1, p.PosEnd(right)); // with one more point on the left
            }
            if (p.PosEnd(right) < p.end())
            {
              areas[3] = simpson_(p.PosBegin(left), p.PosEnd(right) + 1); // with one more point on the right
            }
            UInt valids = 0;
            for (const auto& area : areas)
            {
              if (area != -1.0)
              {
                pa.area += area;
                ++valids;
              }
            }
            pa.area /= valids;
          }
        }
      }
      else if (integration_type_ == INTEGRATION_TYPE_INTENSITYSUM)
      {
        pa.area = compute_peak_area_intensity_sum(left, right);
      }
      else
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Please set a valid value for the parameter \"integration_type\".");
      }
    
      return pa;
    }



    template <typename PeakContainerT>
    PeakBackground estimateBackground_(
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
          for (auto it = p.PosBegin(left); it != p.PosEnd(right); ++it) //OMS_CODING_TEST_EXCLUDE
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

    /**
      @brief Simpson's rule algorithm

      This implementation expects an odd number of points. The formula used supports
      unequally spaced points.

      @note Make sure the container (chromatogram or spectrum) is sorted with respect to position (RT or m/z).

      @warning An odd number of points is expected!

      @param[in] it_begin The iterator to the first point
      @param[in] it_end The iterator to the past-the-last point
      @return The computed area
    */
    template <typename PeakContainerConstIteratorT>
    double simpson_(PeakContainerConstIteratorT it_begin, PeakContainerConstIteratorT it_end) const
    {
      double integral = 0.0;
      for (auto it = it_begin + 1; it < it_end - 1; it = it + 2) //OMS_CODING_TEST_EXCLUDE
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
    PeakShapeMetrics calculatePeakShapeMetrics_(
      const PeakContainerT& pc, double left, double right,
      const double peak_height, const double peak_apex_pos
    ) const
    {
      PeakShapeMetrics psm;

      if (pc.empty()) return psm; // return all '0'

      // enforce order: left <= peakapex <= right
      if (!(left <= peak_apex_pos && peak_apex_pos <= right)) throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);

      PeakContainerT emg_pc;
      const PeakContainerT& p = EMGPreProcess_(pc, emg_pc, left, right);
    
      typename PeakContainerT::ConstIterator it_PosBegin_l = p.PosBegin(left);
      typename PeakContainerT::ConstIterator it_PosEnd_apex = p.PosBegin(peak_apex_pos); // if peak_apex_pos is correct, this will get the underlying iterator
      typename PeakContainerT::ConstIterator it_PosEnd_r = p.PosEnd(right); // past the end. Do not dereference (might be the true .end())
      for (auto it = it_PosBegin_l; it != it_PosEnd_r; ++it) //OMS_CODING_TEST_EXCLUDE
      {
        // points across the peak
        ++(psm.points_across_baseline);
        if (it->getIntensity() >= 0.5 * peak_height)
        {
          ++(psm.points_across_half_height);
        }
      }
      // positions at peak heights
      psm.start_position_at_5 = findPosAtPeakHeightPercent_(it_PosBegin_l, it_PosEnd_apex, p.end(), peak_height, 0.05, true);
      psm.start_position_at_10 = findPosAtPeakHeightPercent_(it_PosBegin_l, it_PosEnd_apex, p.end(), peak_height, 0.1, true);
      psm.start_position_at_50 = findPosAtPeakHeightPercent_(it_PosBegin_l, it_PosEnd_apex, p.end(), peak_height, 0.5, true);
      psm.end_position_at_5 = findPosAtPeakHeightPercent_(it_PosEnd_apex, it_PosEnd_r, p.end(), peak_height, 0.05, false);
      psm.end_position_at_10 = findPosAtPeakHeightPercent_(it_PosEnd_apex, it_PosEnd_r, p.end(), peak_height, 0.1, false);
      psm.end_position_at_50 = findPosAtPeakHeightPercent_(it_PosEnd_apex, it_PosEnd_r, p.end(), peak_height, 0.5, false);
      // peak widths
      psm.width_at_5 = psm.end_position_at_5 - psm.start_position_at_5;
      psm.width_at_10 = psm.end_position_at_10 - psm.start_position_at_10;
      psm.width_at_50 = psm.end_position_at_50 - psm.start_position_at_50;
      psm.total_width = (p.PosEnd(right) - 1)->getPos() - p.PosBegin(left)->getPos();
      psm.slope_of_baseline = (p.PosEnd(right) - 1)->getIntensity() - p.PosBegin(left)->getIntensity();
      if (peak_height != 0.0) // avoid division by zero
      {
        psm.baseline_delta_2_height = psm.slope_of_baseline / peak_height;
      }
      // Source of tailing_factor and asymmetry_factor formulas:
      // USP 40 - NF 35 The United States Pharmacopeia and National Formulary - Supplementary

      // Can only compute if start and peak apex are different 
      if (psm.start_position_at_5 != peak_apex_pos)
      {
        psm.tailing_factor = psm.width_at_5 / (2*(peak_apex_pos - psm.start_position_at_5));
      }
      if (psm.start_position_at_10 != peak_apex_pos)
      {
        psm.asymmetry_factor = (psm.end_position_at_10 - peak_apex_pos) / (peak_apex_pos - psm.start_position_at_10);
      }
      return psm;
    }



    /**
      @brief Find the position (RT/MZ) at a given percentage of peak's height

      @note The method expects that the iterators span half of the peak's width.
      Examples:
      - Left half case: the range would be [leftMostPt, peakApexPos)
      - Right half case: the range would be [peakApexPos + 1, rightMostPt + 1)

      @note The method assumes a convex peak. If 5%, 10%, or 50% peak heights are not found on either side of the peak,
      the closest left (for left peak height percentages) and closest right (for right peak height percentages) will be used.

      @param[in] it_left The iterator to the first point (must not be past the end)
      @param[in] it_right The iterator to the last point (might be past the end)
      @param[in] it_end The end-iterator of the container
      @param[in] peak_height The peak's height
      @param[in] percent At which percentage of the peak height we want to find the position (common values: 0.05, 0.1, 0.5)
      @param[in] is_left_half According to which half of the peak, the algorithm proceeds to the correct direction

      @return The position found
    */
    template <typename PeakContainerConstIteratorT>
    double findPosAtPeakHeightPercent_(
      PeakContainerConstIteratorT it_left,  // must not be past the end
      PeakContainerConstIteratorT it_right, // might be past the end
      PeakContainerConstIteratorT it_end,   // definitely past-the-end
      const double peak_height,
      const double percent,
      const bool is_left_half) const
    {
      // no points in range
      if (it_left == it_end) throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);

      // only one point in range
      if (it_left == it_right) return it_left->getPos();

      const double percent_intensity = peak_height * percent;
      PeakContainerConstIteratorT closest;
      if (is_left_half)
      {
        closest = it_left;
        for (
          PeakContainerConstIteratorT it = it_left;
          it < it_right && it->getIntensity() <= percent_intensity;
          closest = it++
        ) {}
      }
      else // right half; search from right to left
      {
        closest = it_right - 1; // make sure we can deference it
        for (
          PeakContainerConstIteratorT it = it_right - 1;
          it >= it_left && it->getIntensity() <= percent_intensity;
          closest = it--
        ) {}
      }
      return closest->getPos();
    }

private:
    /** @name Parameters
      The user is supposed to select a value for these parameters.
      By default, the integration_type_ is "intensity_sum" and the baseline_type_ is "base_to_base".
    */
    ///@{
    /**
      The integration technique to use in integratePeak() and estimateBackground().
      Possible values are: "trapezoid", "simpson", "intensity_sum".
    */
    String integration_type_ = INTEGRATION_TYPE_INTENSITYSUM;
    /**
      The baseline type to use in estimateBackground().
      Possible values are: "vertical_division_max", "vertical_division_min", "base_to_base".
    */
    String baseline_type_ = BASELINE_TYPE_BASETOBASE;
    ///@}

    /// Enable/disable EMG peak model fitting
    bool fit_EMG_;
    EmgGradientDescent emg_;


    /**
      @brief Fit the peak to the EMG model

      The fitting process happens only if `fit_EMG_` is true. `left` and `right`
      are updated accordingly.

      @tparam PeakContainerT Either a MSChromatogram or a MSSpectrum
      @param[in] pc Input peak
      @param[out] emg_pc Will possibly contain the processed peak
      @param[in] left RT or MZ value of the first point of interest
      @param[in] right RT or MZ value of the first point of interest
      @return A const reference to `emg_pc` if the fitting is executed, `pc` otherwise.
    */
    template <typename PeakContainerT>
    const PeakContainerT& EMGPreProcess_(
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
  };
}
