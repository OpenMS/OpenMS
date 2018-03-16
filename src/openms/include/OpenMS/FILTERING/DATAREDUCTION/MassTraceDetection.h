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
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_DATAREDUCTION_MASSTRACEDETECTION_H
#define OPENMS_FILTERING_DATAREDUCTION_MASSTRACEDETECTION_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{

  /**
    @brief A mass trace extraction method that gathers peaks similar in m/z and moving along retention time.

    Peaks of a @ref MSExperiment are sorted by their intensity and stored in a
    list of potential chromatographic apex positions. Only peaks that are above
    the noise threshold (user-defined) are analyzed and only peaks that are n
    times above this minimal threshold are considered as apices. This saves
    computational resources and decreases the noise in the resulting output.
    
    Starting with these, mass traces are extended in- and decreasingly in
    retention time. During this extension phase, the centroid m/z is computed
    on-line as an intensity-weighted mean of peaks.
    
    The extension phase ends when either the frequency of gathered peaks drops
    below a threshold (min_sample_rate, see @ref MassTraceDetection parameters)
    or when the number of missed scans exceeds a threshold
    (trace_termination_outliers, see @ref MassTraceDetection parameters).

    Finally, only mass traces that pass a filter (a certain minimal and maximal
    length as well as having the minimal sample rate criterion fulfilled) get
    added to the result.

    @htmlinclude OpenMS_MassTraceDetection.parameters

    @ingroup Quantitation
  */
  class OPENMS_DLLAPI MassTraceDetection :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Default constructor
    MassTraceDetection();

    /// Default destructor
    ~MassTraceDetection() override;

    /** @name Helper methods
    */

    /// Allows the iterative computation of the intensity-weighted mean of a mass trace's centroid m/z.
    void updateIterativeWeightedMeanMZ(const double &, const double &, double &, double &, double &);

    /** @name Main computation methods
    */

    /// Main method of MassTraceDetection. Extracts mass traces of a @ref MSExperiment and gathers them into a vector container.
    void run(const PeakMap &, std::vector<MassTrace> &);

    /// Invokes the run method (see above) on merely a subregion of a @ref MSExperiment map.
    void run(PeakMap::ConstAreaIterator & begin, PeakMap::ConstAreaIterator & end, std::vector<MassTrace> & found_masstraces);

    /** @name Private methods and members 
    */
protected:
    void updateMembers_() override;

private:

    typedef std::multimap<double, std::pair<Size, Size> > MapIdxSortedByInt;

    /// The internal run method
    void run_(const MapIdxSortedByInt& chrom_apices,
              const Size peak_count, 
              const PeakMap & work_exp,
              const std::vector<Size>& spec_offsets,
              std::vector<MassTrace> & found_masstraces);

    // parameter stuff
    double mass_error_ppm_;
    double noise_threshold_int_;
    double chrom_peak_snr_;
    MassTrace::MT_QUANTMETHOD quant_method_;

    String trace_termination_criterion_;
    Size trace_termination_outliers_;
    double min_sample_rate_;
    double min_trace_length_;
    double max_trace_length_;

    bool reestimate_mt_sd_;
  };
}

#endif // OPENMS_FILTERING_DATAREDUCTION_MASSTRACEDETECTION_H
