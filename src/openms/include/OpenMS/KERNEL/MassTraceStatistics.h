// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/OpenMSConfig.h>

#include <span>

namespace OpenMS
{
  /**
    @brief Stateless statistics over a range of mass trace peaks.

    The functions take a contiguous range of Peak2D sorted by retention time
    (e.g. a MassTrace, a std::vector<Peak2D>, or a part of one) and return
    their result instead of storing it. This allows evaluating a candidate
    (sub-)trace without building a MassTrace for it first, e.g.

    @code
    std::span<const Peak2D> peaks(trace.getPeaks());
    double mz = MassTraceStatistics::weightedMeanMZ(peaks.subspan(first, count));
    @endcode

    MassTrace uses these functions for its update*(), estimateFWHM() and
    compute*() members and caches the results.

    Functions taking @p smoothed expect one smoothed intensity per peak, i.e.
    <tt>smoothed.size() == peaks.size()</tt>.

    @ingroup Kernel
  */
  namespace MassTraceStatistics
  {
    /// A contiguous range of peaks, sorted by retention time
    using PeakRange = std::span<const Peak2D>;

    /// One (e.g. smoothed) intensity per peak of a PeakRange
    using IntensityRange = std::span<const double>;

    /// Result of estimateFWHM()
    struct FWHM
    {
      double fwhm = 0.0;  ///< full width at half maximum (in RT units); 0 if the apex is the first or last peak
      Size start_idx = 0; ///< index of the first peak of the FWHM region (inclusive)
      Size end_idx = 0;   ///< index of the last peak of the FWHM region (inclusive)
    };

    /** @name Centroid m/z
    */
    ///@{

    /// Mean of the m/z values.
    /// @throw Exception::InvalidValue if @p peaks is empty
    OPENMS_DLLAPI double meanMZ(PeakRange peaks);

    /// Median of the m/z values.
    /// @throw Exception::InvalidValue if @p peaks is empty
    OPENMS_DLLAPI double medianMZ(PeakRange peaks);

    /// Intensity-weighted mean of the m/z values.
    /// @throw Exception::InvalidValue if @p peaks is empty, or holds several peaks whose intensities sum up to zero
    OPENMS_DLLAPI double weightedMeanMZ(PeakRange peaks);

    /// Intensity-weighted standard deviation of the m/z values around @p centroid_mz.
    /// @throw Exception::InvalidValue if @p peaks is empty or the intensities sum up to zero
    OPENMS_DLLAPI double weightedMZsd(PeakRange peaks, double centroid_mz);
    ///@}

    /** @name Centroid RT
    */
    ///@{

    /// Intensity-weighted mean of the RTs (trapezoidal rule; plain mean if the weights carry no information).
    /// @throw Exception::InvalidValue if @p peaks is empty
    OPENMS_DLLAPI double weightedMeanRT(PeakRange peaks);

    /// Mean of the RTs, weighted by the positive @p smoothed intensities.
    /// @throw Exception::InvalidValue if @p smoothed is empty, or @p peaks holds several peaks and no smoothed intensity is positive
    OPENMS_DLLAPI double smoothedWeightedMeanRT(PeakRange peaks, IntensityRange smoothed);

    /// RT of the peak with the highest @p smoothed intensity.
    /// @throw Exception::InvalidValue if @p smoothed is empty, or @p peaks holds several peaks and no smoothed intensity is positive
    OPENMS_DLLAPI double smoothedMaxRT(PeakRange peaks, IntensityRange smoothed);

    /// Median of the RTs.
    /// @throw Exception::InvalidValue if @p peaks is empty
    OPENMS_DLLAPI double medianRT(PeakRange peaks);
    ///@}

    /** @name Apex and peak width
    */
    ///@{

    /// Index of the most intense peak (the first one, if several share the maximum).
    /// @throw Exception::InvalidValue if @p peaks is empty
    OPENMS_DLLAPI Size apexIndex(PeakRange peaks);

    /// Index of the peak with the highest @p smoothed intensity (the first one, if several share the maximum).
    /// @throw Exception::InvalidValue if @p smoothed or @p peaks is empty
    OPENMS_DLLAPI Size apexIndex(PeakRange peaks, IntensityRange smoothed);

    /**
      @brief Estimates the FWHM of the chromatographic peak around the apex.

      The RTs at half maximum are linearly interpolated between the two peaks
      enclosing the half maximum on either flank. If a flank does not drop
      below the half maximum, the RT of its outermost peak is used.

      @pre @p peaks is sorted by RT
      @throw Exception::InvalidValue if @p peaks is empty
    */
    OPENMS_DLLAPI FWHM estimateFWHM(PeakRange peaks);

    /// As estimateFWHM(PeakRange), but based on the @p smoothed intensities.
    /// @throw Exception::InvalidValue if @p smoothed or @p peaks is empty
    OPENMS_DLLAPI FWHM estimateFWHM(PeakRange peaks, IntensityRange smoothed);
    ///@}

    /** @name Intensities and areas
    */
    ///@{

    /// Sum of the peak intensities
    OPENMS_DLLAPI double intensitySum(PeakRange peaks);

    /// Highest peak intensity (0 for an empty range)
    OPENMS_DLLAPI double maxIntensity(PeakRange peaks);

    /// Highest of the @p intensities (0 for an empty range)
    OPENMS_DLLAPI double maxIntensity(IntensityRange intensities);

    /// Median of the peak intensities.
    /// @pre @p peaks is not empty
    OPENMS_DLLAPI double medianIntensity(PeakRange peaks);

    /// Area under the peaks over RT (trapezoidal rule); 0 for less than two peaks
    OPENMS_DLLAPI double peakArea(PeakRange peaks);

    /// Area under the @p smoothed intensities over the RTs of @p peaks (trapezoidal rule); 0 for less than two peaks
    OPENMS_DLLAPI double peakArea(PeakRange peaks, IntensityRange smoothed);
    ///@}
  } // namespace MassTraceStatistics
} // namespace OpenMS
