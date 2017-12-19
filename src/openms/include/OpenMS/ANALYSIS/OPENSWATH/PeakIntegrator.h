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
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  /**
    @brief Compute the area, background and shape metrics of a peak.

    The containers supported by the methods are MSChromatogram and MSSpectrum.
  */
  class OPENMS_DLLAPI PeakIntegrator :
    public DefaultParamHandler
  {
public:
    /// Constructor
    PeakIntegrator();
    /// Destructor
    virtual ~PeakIntegrator();

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

    /**
      @brief Compute the area of a peak contained in a MSChromatogram.

      The value of integration_type_ decides which integration technique to use:
      - "trapezoid" for the trapezoidal rule
      - "simpson" for the Simpson's rule (for unequally spaced points, Shklov, 1960)
      - "intensity_sum" for the simple sum of the intensities

      @note Make sure the chromatogram is sorted with respect to retention time.

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
    void updateMembers_();

    template <typename PeakContainerT>
    PeakArea integratePeak_(
      const PeakContainerT& p, const double left, const double right
    ) const;

    template <typename PeakContainerT>
    PeakBackground estimateBackground_(
      const PeakContainerT& p, const double left, const double right,
      const double peak_apex_pos
    ) const;

    template <typename PeakContainerConstIteratorT>
    double simpson_(PeakContainerConstIteratorT it_begin, PeakContainerConstIteratorT it_end) const;

    template <typename PeakContainerT>
    PeakShapeMetrics calculatePeakShapeMetrics_(
      const PeakContainerT& p, const double left, const double right,
      const double peak_height, const double peak_apex_pos
    ) const;

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
      Possible values are: "vertical_division", "base_to_base".
    */
    String baseline_type_ = BASELINE_TYPE_BASETOBASE;
    ///@}

    /** @name Helper methods
      The Simpson's rule implementations for an odd number of unequally spaced points.
    */
    ///@{
    /**
      @brief Simpson's rule algorithm

      This implementation expects an odd number of points. The formula used supports
      unequally spaced points.

      @note Make sure the chromatogram is sorted with respect to retention time.

      @warning An odd number of points is expected!

      @param[in] it_begin The iterator to the first point
      @param[in] it_end The iterator to the past-the-last point
      @return The computed area
    */
    double simpson(MSChromatogram::ConstIterator it_begin, MSChromatogram::ConstIterator it_end) const;
    /**
      @brief Simpson's rule algorithm

      This implementation expects an odd number of points. The formula used supports
      unequally spaced points.

      @note Make sure the spectrum is sorted with respect to mass-to-charge ratio.

      @warning An odd number of points is expected!

      @param[in] it_begin The iterator to the first point
      @param[in] it_end The iterator to the past-the-last point
      @return The computed area
    */
    double simpson(MSSpectrum::ConstIterator it_begin, MSSpectrum::ConstIterator it_end) const;
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
    /// Baseline type: vertical division
    static constexpr const char* BASELINE_TYPE_VERTICALDIVISION = "vertical_division";
    ///@}
  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_PEAKINTEGRATOR_H
