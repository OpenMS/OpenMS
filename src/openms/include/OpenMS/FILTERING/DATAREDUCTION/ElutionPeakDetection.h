// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------


#ifndef OPENMS_FILTERING_DATAREDUCTION_ELUTIONPEAKDETECTION_H
#define OPENMS_FILTERING_DATAREDUCTION_ELUTIONPEAKDETECTION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MassTrace.h>

namespace OpenMS
{
  /**
    @brief Extracts chromatographic peaks from a mass trace.

    Mass traces may consist of several consecutively (partly overlapping)
    eluting peaks, e.g., stemming from isomeric compounds with exactly the
    same mass but different retentional behaviour. This method first applies
    smoothing on the mass trace's intensities, then detects local
    minima/maxima in order to separate the chromatographic peaks from each
    other. This results in a vector that gathers the split mass traces
    (see @ref ElutionPeakDetection parameters).

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
    virtual ~ElutionPeakDetection();

    /** @brief Extracts chromatographic peaks from a single MassTrace and stores the splits into a vector of new mass traces.
     *
     * @note Smoothed intensities are added to @p mt_vec
     *
     * @param mt Input mass trace
     * @param single_mtraces Output single mass traces (detected peaks)
     *
    */
    void detectPeaks(MassTrace& mt, std::vector<MassTrace>& single_mtraces);

    /** @brief Extracts chromatographic peaks from multiple MassTraces and stores the splits into a vector of new mass traces.
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
     * This function computes local extrema on a given input mass trace
     *
     * @param tr Input mass trace
     * @param num_neighboring_peaks How many data points are expected to belong to a peak (used to split traces and find maxima)
     * @param chrom_maxes Output of maxima (gets cleared)
     * @param chrom_mins Output of minima (gets cleared)
     *
    */
    void findLocalExtrema(const MassTrace& tr, const Size& num_neighboring_peaks,
                          std::vector<Size>& chrom_maxes, std::vector<Size>& chrom_mins);

    /// adds smoothed_intensities to internal data of @p mt
    void smoothData(MassTrace& mt, int win_size) const;

protected:
    virtual void updateMembers_();

private:
    double chrom_fwhm_;
    // double scan_time_;
    double chrom_peak_snr_;
    //double noise_threshold_int_;
    //double sample_rate_;

    double min_fwhm_;
    double max_fwhm_;
    // double min_trace_length_;
    // double max_trace_length_;

    String pw_filtering_;
    bool mt_snr_filtering_;

    void detectElutionPeaks_(MassTrace&, std::vector<MassTrace>&);
  };

} // namespace OpenMS

#endif // OPENMS_FILTERING_DATAREDUCTION_ELUTIONPEAKDETECTION_H
