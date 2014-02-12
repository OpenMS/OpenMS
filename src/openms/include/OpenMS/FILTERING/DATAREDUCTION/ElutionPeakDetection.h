// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

    Mass traces may consist of several consecutively (partly overlapping) eluting peaks, e.g., stemming
    from isomeric compounds with exactly the same mass but different retentional behaviour. This method
    first applies LOWESS smoothing on the mass trace's intensities, then detects local minima/maxima in
    order to separate the chromatographic peaks from each other. This results in a vector that gathers
    the splitted mass traces (see @ref ElutionPeakDetection parameters).

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

    /// Getter & setter methods

    /** @name Main computation methods
        */
    /// Extracts chromatographic peaks from a single MassTrace and stores the splits into a vector of new mass traces.
    void detectPeaks(MassTrace &, std::vector<MassTrace> &);

    /// Applies the aforementioned detection method on a series of mass traces as input.
    void detectPeaks(std::vector<MassTrace> &, std::vector<MassTrace> &);

    void filterByPeakWidth(std::vector<MassTrace> &, std::vector<MassTrace> &);
    DoubleReal computeMassTraceNoise(const MassTrace&);
    DoubleReal computeMassTraceSNR(const MassTrace&);
    DoubleReal computeApexSNR(const MassTrace&);
    void findLocalExtrema(const MassTrace&, const Size&, std::vector<Size>&, std::vector<Size>&);

protected:
    virtual void updateMembers_();

private:
    DoubleReal chrom_fwhm_;
    // DoubleReal scan_time_;
    DoubleReal chrom_peak_snr_;
    DoubleReal noise_threshold_int_;
    DoubleReal sample_rate_;

    DoubleReal min_fwhm_;
    DoubleReal max_fwhm_;
    // DoubleReal min_trace_length_;
    // DoubleReal max_trace_length_;

    String pw_filtering_;
    bool mt_snr_filtering_;

    void detectElutionPeaks_(MassTrace &, std::vector<MassTrace> &);
  };
} // namespace OpenMS
#endif // OPENMS_FILTERING_DATAREDUCTION_ELUTIONPEAKDETECTION_H
