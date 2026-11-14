// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MassTrace.h>

namespace OpenMS
{
  /**
    @brief Extracts chromatographic peaks from a mass trace.

    Mass traces may consist of several consecutively (partly overlapping)
    eluting peaks, e.g., stemming from (almost) isobaric compounds that are
    separated by retention time. Especially in metabolomics, isomeric compounds
    with exactly the same mass but different retentional behaviour may still be
    contained in the same mass trace.
    
    This method first applies smoothing on the mass trace's intensities, then
    detects local minima/maxima in order to separate the chromatographic peaks
    from each other. Detection of maxima is performed on the smoothed
    intensities and uses a fixed peak width (given as parameter chrom_fwhm)
    within which only a single maximum is expected. Currently smoothing is
    done using SavitzkyGolay smoothing with a second order polynomial and a
    frame length of the fixed peak width (see smoothData() for how the ends
    of a mass trace are treated).

    Depending on the "width_filtering" parameters, mass traces are filtered by
    length in seconds ("fixed" filter) or by quantile.

    The output of the algorithm is a set of chromatographic
    peaks for each mass trace, i.e. a vector of split mass traces (see @ref
    ElutionPeakDetection parameters).

    In general, a user would want to call the "detectPeaks" functions,
    potentially followed by the "filterByPeakWidth" function.

    @htmlinclude OpenMS_ElutionPeakDetection.parameters

    @ingroup Quantitation
  */
  class OPENMS_DLLAPI ElutionPeakDetection :
    public DefaultParamHandler, public ProgressLogger
  {
public:
    /// Default Constructor
    ElutionPeakDetection();

    /// Destructor
    ~ElutionPeakDetection() override;

    /** @brief Extracts chromatographic peaks from a single MassTrace and
     *         stores the resulting split traces in a vector of new mass traces.
     *
     * @note Smoothed intensities are added to @p mt_vec
     *
     * @param[in,out] mt Input mass trace
     * @param[out] single_mtraces Output single mass traces (detected peaks)
     *
    */
    void detectPeaks(MassTrace& mt, std::vector<MassTrace>& single_mtraces);

    /** @brief Extracts chromatographic peaks from multiple MassTraces and
     *         stores the resulting split traces in a vector of new mass traces.
     *
     * @note Smoothed intensities are added to @p mt_vec
     *
     * @param[in,out] mt_vec Input mass traces
     * @param[out] single_mtraces Output single mass traces (detected peaks)
     *
    */
    void detectPeaks(std::vector<MassTrace>& mt_vec, std::vector<MassTrace>& single_mtraces);

    /// Filter out mass traces below lower 5 % quartile and above upper 95 % quartile
    void filterByPeakWidth(std::vector<MassTrace>&, std::vector<MassTrace>&);

    /// Compute noise level (as RMSE of the actual signal and the smoothed signal)
    double computeMassTraceNoise(const MassTrace&);

    /// Compute the signal to noise ratio (estimated by computeMassTraceNoise)
    double computeMassTraceSNR(const MassTrace&);

    /// Compute the signal to noise ratio at the apex (estimated by computeMassTraceNoise)
    double computeApexSNR(const MassTrace&);

    /** @brief Computes local extrema on a mass trace
     *
     * This function computes local extrema on a given input mass trace. It
     * works on the smoothed intensities which must be available at this step.
     * Initially it identifies potential maxima as peaks that have maximum
     * intensity within a range of peak +/- num_neighboring_peaks.
     * All such maxima in the smoothed data get added to the list of maxima.
     * Minima are found through bisection between the maxima.
     *
     * @param[in] tr Input mass trace
     * @param[in] num_neighboring_peaks How many data points are expected to belong
     *                              to a peak, i.e. the expected peak width
     *                              (this is used to split traces and find
     *                              maxima)
     * @param[out] chrom_maxes Output of maxima (gets cleared)
     * @param[out] chrom_mins Output of minima (gets cleared)
     *
     * Returns a vector of indices where a maxima may occur (chrom_maxes) and a
     * vector of indices where a minima may occur (chrom_mins).
     *
     * @note this expects that the input mass trace has been smoothed before
     * and the smoothed intensities are available through tr.getSmoothedIntensities().
     *
    */
    void findLocalExtrema(const MassTrace& tr, const Size& num_neighboring_peaks,
                          std::vector<Size>& chrom_maxes, std::vector<Size>& chrom_mins) const;

    /** @brief Adds smoothed intensities to the internal data of @p mt
     *
     * Savitzky-Golay smoothing with a second order polynomial. The frame is
     * centred on each data point. Towards the ends of the trace, where a
     * full frame does not fit, it narrows to the largest one that does,
     * instead of the fit of the outermost full frame being extrapolated: a
     * mass trace often starts or ends on a steep flank, and an extrapolated
     * fit can put the smoothed maximum on the first or last point of the
     * trace. The second and the last but one point are smoothed with the
     * mean of three points (a three point frame is fitted exactly by a
     * second order polynomial), the first and the last point keep their
     * intensity. Each smoothed intensity is limited to the range of the raw
     * intensities within its frame and to >= 0, so that the fit cannot ring
     * below the lowest (in particular: to zero) or above the highest of them.
     *
     * @param[in,out] mt Mass trace to smooth
     * @param[in] win_size Frame length in data points, the number of points
     *                     one expected peak width spans. Even values are
     *                     incremented by one. A frame of three points (or
     *                     less) does not smooth, and neither does a frame
     *                     that is longer than the trace.
    */
    void smoothData(MassTrace& mt, int win_size) const;

protected:
    void updateMembers_() override;

private:
    double chrom_fwhm_;
    double chrom_peak_snr_;

    double min_fwhm_;
    double max_fwhm_;

    /// Type of width filtering
    std::string pw_filtering_;

    /// Whether to apply S/N filtering
    bool mt_snr_filtering_;

    /// Main function to do the work
    void detectElutionPeaks_(MassTrace&, std::vector<MassTrace>&);
  };

} // namespace OpenMS

