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
    frame length of the fixed peak width.

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
     * @param mt Input mass trace
     * @param single_mtraces Output single mass traces (detected peaks)
     *
    */
    void detectPeaks(MassTrace& mt, std::vector<MassTrace>& single_mtraces);

    /** @brief Extracts chromatographic peaks from multiple MassTraces and
     *         stores the resulting split traces in a vector of new mass traces.
     *
     * @note Smoothed intensities are added to @p mt_vec
     *
     * @param mt_vec Input mass traces
     * @param single_mtraces Output single mass traces (detected peaks)
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
     * @param tr Input mass trace
     * @param num_neighboring_peaks How many data points are expected to belong
     *                              to a peak, i.e. the expected peak width
     *                              (this is used to split traces and find
     *                              maxima)
     * @param chrom_maxes Output of maxima (gets cleared)
     * @param chrom_mins Output of minima (gets cleared)
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

    /// adds smoothed_intensities to internal data of @p mt
    void smoothData(MassTrace& mt, int win_size) const;

protected:
    void updateMembers_() override;

private:
    double chrom_fwhm_;
    double chrom_peak_snr_;

    double min_fwhm_;
    double max_fwhm_;

    /// Type of width filtering
    String pw_filtering_;

    /// Whether to apply S/N filtering
    bool mt_snr_filtering_;

    /// Main function to do the work
    void detectElutionPeaks_(MassTrace&, std::vector<MassTrace>&);
  };

} // namespace OpenMS

