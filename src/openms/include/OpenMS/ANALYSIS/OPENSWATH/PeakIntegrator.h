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

#pragma once

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

    /// To test private and protected methods
    friend class PeakIntegrator_friend;

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

    /**
      @brief Fit the given peak (a MSChromatogram) to the EMG peak model

      The method is able to recapitulate the actual peak area of saturated or cutoff peaks.
      In addition, the method is able to fine tune the peak area of well acquired peaks.
      The output is a reconstruction of the input peak. Additional points are often added
      to produce a peak with similar intensities on boundaries' points.

      Metadata will be added to the output peak, containing the optimal parameters
      for the EMG peak model. This information will be found in a `FloatDataArray`
      of name "emg_parameters", with the parameters being saved in the following
      order (from index 0 to 3): amplitude `h`, mean `mu`, standard deviation `sigma`,
      exponent relaxation time `tau`.

      @note All optimal gradient descent parameters are currently hard coded to allow for a simplified user interface

      @note Cutoff peak: The intensities of the left and right baselines are not equal

      @note Saturated peak: The maximum intensity of the peak is lower than expected due to saturation of the detector

      Inspired by the results found in:
      Yuri Kalambet, Yuri Kozmin, Ksenia Mikhailova, Igor Nagaev, Pavel Tikhonov
      Reconstruction of chromatographic peaks using the exponentially modified Gaussian function

      @param[in] input_peak Input peak
      @param[out] output_peak Output peak
    */
    void fitEMGPeakModel(
      const MSChromatogram& input_peak,
      MSChromatogram& output_peak
    ) const;

    /**
      @brief Fit the given peak (a MSSpectrum) to the EMG peak model

      The method is able to recapitulate the actual peak area of saturated or cutoff peaks.
      In addition, the method is able to fine tune the peak area of well acquired peaks.
      The output is a reconstruction of the input peak. Additional points are often added
      to produce a peak with similar intensities on boundaries' points.

      Metadata will be added to the output peak, containing the optimal parameters
      for the EMG peak model. This information will be found in a `FloatDataArray`
      of name "emg_parameters", with the parameters being saved in the following
      order (from index 0 to 3): amplitude `h`, mean `mu`, standard deviation `sigma`,
      exponent relaxation time `tau`.

      @note All optimal gradient descent parameters are currently hard coded to allow for a simplified user interface

      @note Cutoff peak: The intensities of the left and right baselines are not equal

      @note Saturated peak: The maximum intensity of the peak is lower than expected due to saturation of the detector

      Inspired by the results found in:
      Yuri Kalambet, Yuri Kozmin, Ksenia Mikhailova, Igor Nagaev, Pavel Tikhonov
      Reconstruction of chromatographic peaks using the exponentially modified Gaussian function

      @param[in] input_peak Input peak
      @param[out] output_peak Output peak
    */
    void fitEMGPeakModel(
      const MSSpectrum& input_peak,
      MSSpectrum& output_peak
    ) const;

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

    /**
      @brief Find the position (RT/MZ) at a given percentage of peak's height

      @note The method expects that the iterators span half of the peak's width.
      Examples:
      - Left half case: the range would be [leftMostPt, peakApexPos)
      - Right half case: the range would be [peakApexPos + 1, rightMostPt + 1)

      @note The method assumes a convex peak. If 5%, 10%, or 50% peak heights are not found on either side of the peak,
      the closest left (for left peak height percentages) and closest right (for right peak height percentages) will be used.

      @param[in] it_begin The iterator to the first point
      @param[in] it_end The iterator to past-the-last point
      @param[in] peak_height The peak's height
      @param[in] percent At which percentage of the peak height we want to find the position (common values: 0.05, 0.1, 0.5)
      @param[in] is_left_half According to which half of the peak, the algorithm proceeds to the correct direction

      @return The position found
    */
    template <typename PeakContainerConstIteratorT>
    double findPosAtPeakHeightPercent_(
      PeakContainerConstIteratorT it_begin,
      PeakContainerConstIteratorT it_end,
      const double peak_height,
      const double percent,
      const bool is_left_half
    ) const;

    /**
      @brief Given a peak, extract a training set to be used with the gradient descent algorithm

      The algorithm tries to select only those points that can help in finding the optimal
      parameters with gradient descent. The decision of which points to skip is based on the
      derivatives between consecutive points.

      It first selects all those points whose intensity is below a certain value (`intensity_threshold`).
      Then, the derivatives of all the remaining points are computed. Based on the results,
      the algorithm selects those points that present a high enough derivative.
      Once a low value is found, the algorithm stops taking points from that side.
      It then repeats the same procedure on the other side of the peak.
      The goal is to limit the inclusion of saturated or spurious points near the
      peak apex during training.

      @throw Exception::SizeUnderflow if the input has less than 2 elements

      @param[in] xs Positions
      @param[in] ys Intensities
      @param[out] TrX Extracted training set positions
      @param[out] TrY Extracted training set intensities
    */
    void extractTrainingSet(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      std::vector<double>& TrX,
      std::vector<double>& TrY
    ) const;

    /**
      @brief Compute the boundary for the mean (`mu`) parameter in gradient descent

      Together with the value returned by computeInitialMean(), this method
      decides the minimum and maximum value that `mu` can assume during iterations
      of the gradient descent algorithm.
      The value is based on the width of the peak.

      @param[in] xs Positions

      @return The maximum distance from the precomputed initial mean in the gradient descent algorithm
    */
    double computeMuMaxDistance(const std::vector<double>& xs) const;

    /**
      @brief Compute an estimation of the mean of a peak

      The method computes the middle point on different levels of intensity of the peak.
      The returned mean is the average of these middle points.

      @throw Exception::SizeUnderflow if the input is empty

      @param[in] xs Positions
      @param[in] ys Intensities

      @return The peak's estimated mean
    */
    double computeInitialMean(
      const std::vector<double>& xs,
      const std::vector<double>& ys
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
      Possible values are: "vertical_division_max", "vertical_division_min", "base_to_base".
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

    /**
      @brief The implementation of the gradient descent algorithm for the EMG peak model

      @param[in] xs Positions
      @param[in] ys Intensities
      @param[out] best_h `h` (amplitude) parameter
      @param[out] best_mu `mu` (mean) parameter
      @param[out] best_sigma `sigma` (standard deviation) parameter
      @param[out] best_tau `tau` (exponent relaxation time) parameter

      @return The number of iterations necessary to reach the best values for the parameters
    */
    UInt emg_gradient_descent(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      double& best_h,
      double& best_mu,
      double& best_sigma,
      double& best_tau
    ) const;

    /**
      @brief Apply the iRprop+ algorithm for gradient descent

      Reference:
      Christian Igel and Michael Hüsken. Improving the Rprop Learning Algorithm.
      Second International Symposium on Neural Computation (NC 2000), pp. 115-121, ICSC Academic Press, 2000

      @param[in] prev_diff_E_param The cost of the partial derivative of E with
      respect to the given parameter, at the previous iteration of gradient descent
      @param[in,out] diff_E_param The cost of the partial derivative of E with
      respect to the given parameter, at the current iteration
      @param[in,out] param_lr The learning rate for the given parameter
      @param[in,out] param_update The amount to add/remove to/from `param`
      @param[in,out] param The parameter for which the algorithm tries speeding the convergence to a minimum
      @param[in] current_E The current cost E
      @param[in] previous_E The previous cost E
    */
    void iRpropPlus(
      const double prev_diff_E_param,
      double& diff_E_param,
      double& param_lr,
      double& param_update,
      double& param,
      const double current_E,
      const double previous_E
    ) const;

    /**
      @brief Compute the cost given by loss function E

      Needed by the gradient descent algorithm.
      The mean squared error is used as the loss function E.

      @param[in] xs Positions
      @param[in] ys Intensities
      @param[in] h Amplitude
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time

      @return The computed cost
    */
    double Loss_function(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const;

    /**
      @brief Compute the cost given by the partial derivative of the loss function E,
      with respect to `h` (the amplitude)

      Needed by the gradient descent algorithm.

      @param[in] xs Positions
      @param[in] ys Intensities
      @param[in] h Amplitude
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time

      @return The computed cost
    */
    double E_wrt_h(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const;

    /**
      @brief Compute the cost given by the partial derivative of the loss function E,
      with respect to `mu` (the mean)

      Needed by the gradient descent algorithm.

      @param[in] xs Positions
      @param[in] ys Intensities
      @param[in] h Amplitude
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time

      @return The computed cost
    */
    double E_wrt_mu(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const;

    /**
      @brief Compute the cost given by the partial derivative of the loss function E,
      with respect to `sigma` (the standard deviation)

      Needed by the gradient descent algorithm.

      @param[in] xs Positions
      @param[in] ys Intensities
      @param[in] h Amplitude
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time

      @return The computed cost
    */
    double E_wrt_sigma(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const;

    /**
      @brief Compute the cost given by the partial derivative of the loss function E,
      with respect to `tau` (the exponent relaxation time)

      Needed by the gradient descent algorithm.

      @param[in] xs Positions
      @param[in] ys Intensities
      @param[in] h Amplitude
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time

      @return The computed cost
    */
    double E_wrt_tau(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const;

    /**
      @brief Compute EMG's z parameter

      The value of z decides which formula is to be used during EMG function computation.
      Z values in the following ranges will each use a different EMG formula to
      avoid numerical instability and potential numerical overflow:
      (-inf, 0), [0, 6.71e7], (6.71e7, +inf)

      Reference:
      Kalambet, Y.; Kozmin, Y.; Mikhailova, K.; Nagaev, I.; Tikhonov, P. (2011).
      "Reconstruction of chromatographic peaks using the exponentially modified
      Gaussian function". Journal of Chemometrics. 25 (7): 352.

      @param[in] x Position
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time

      @return The computed parameter z
    */
    double compute_z(
      const double x,
      const double mu,
      const double sigma,
      const double tau
    ) const;

    /**
      @brief Compute the EMG function on a set of points

      If class parameter `compute_additional_points` is `"true"`, the algorithm
      will detect which side of the peak is cutoff and add points to it.

      @param[in] xs Positions
      @param[in] h Amplitude
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time
      @param[out] out_xs The output positions
      @param[out] out_ys The output intensities
    */
    void emg_vector(
      const std::vector<double>& xs,
      const double h,
      const double mu,
      const double sigma,
      const double tau,
      std::vector<double>& out_xs,
      std::vector<double>& out_ys
    ) const;

    /**
      @brief Compute the EMG function on a single point

      @param[in] x Position
      @param[in] h Amplitude
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time

      @return The estimated intensity for the given input point
    */
    double emg_point(
      const double x,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const;

    template <typename PeakContainerT>
    void fitEMGPeakModel_(
      const PeakContainerT& input_peak,
      PeakContainerT& output_peak
    ) const;

    /// Alias for OpenMS::Constants:PI
    const double PI = OpenMS::Constants::PI;

    /**
      Level of debug information to print to the terminal
      Valid values are: 0, 1, 2
      Higher values mean more information
    */
    UInt print_debug_;

    /// Maximum number of gradient descent iterations in `fitEMGPeakModel()`
    UInt max_gd_iter_;

    /**
      Whether additional points should be added when fitting EMG peak model,
      particularly useful with cutoff peaks
    */
    bool compute_additional_points_;
  };

  class PeakIntegrator_friend
  {
public:
    PeakIntegrator_friend() = default;
    ~PeakIntegrator_friend() = default;

    double Loss_function(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const
    {
      return peakIntegrator.Loss_function(xs, ys, h, mu, sigma, tau);
    }

    double computeMuMaxDistance(const std::vector<double>& xs) const
    {
      return peakIntegrator.computeMuMaxDistance(xs);
    }

    void extractTrainingSet(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      std::vector<double>& TrX,
      std::vector<double>& TrY
    ) const
    {
      peakIntegrator.extractTrainingSet(xs, ys, TrX, TrY);
    }

    double computeInitialMean(
      const std::vector<double>& xs,
      const std::vector<double>& ys
    ) const
    {
      return peakIntegrator.computeInitialMean(xs, ys);
    }

    void iRpropPlus(
      const double prev_diff_E_param,
      double& diff_E_param,
      double& param_lr,
      double& param_update,
      double& param,
      const double current_E,
      const double previous_E
    ) const
    {
      peakIntegrator.iRpropPlus(
        prev_diff_E_param, diff_E_param, param_lr,
        param_update, param, current_E, previous_E
      );
    }

    double compute_z(
      const double x,
      const double mu,
      const double sigma,
      const double tau
    ) const
    {
      return peakIntegrator.compute_z(x, mu, sigma, tau);
    }

    void emg_vector(
      const std::vector<double>& xs,
      const double h,
      const double mu,
      const double sigma,
      const double tau,
      std::vector<double>& out_xs,
      std::vector<double>& out_ys
    ) const
    {
      peakIntegrator.emg_vector(xs, h, mu, sigma, tau, out_xs, out_ys);
    }

    double emg_point(
      const double x,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const
    {
      return peakIntegrator.emg_point(x, h, mu, sigma, tau);
    }

    PeakIntegrator peakIntegrator;
  };
}

