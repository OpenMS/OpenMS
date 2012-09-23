// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Erhan Kenar $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_DATAREDUCTION_MASSTRACEDETECTION_H
#define OPENMS_FILTERING_DATAREDUCTION_MASSTRACEDETECTION_H

#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{
/**
  @brief A mass trace extraction method that gathers peaks similar in m/z and moving along retention time.

  Peaks of a @ref MSExperiment are sorted by their intensity and stored in a list of potential chromatographic
  apex positions. Starting with these, mass traces are extended in- and decreasingly in retention
  time. During this extension phase, the centroid m/z is computed on-line as an intensity-weighted mean of
  peaks. The extension phase ends when the frequency of gathered peaks drops below a
  threshold (min_sample_rate, see @ref MassTraceDetection parameters).

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
    virtual ~MassTraceDetection();

    /** @name Helper methods
        */
    /// Allows the iterative computation of the intensity-weighted mean of a mass trace's centroid m/z.
    void updateIterativeWeightedMeanMZ(const DoubleReal &, const DoubleReal &, DoubleReal &, DoubleReal &, DoubleReal &);

    /// Computes a rough estimate of the average peak width of the experiment (median) and an estimate of a lower and upper bound for the peak width (+/-2*MAD, median of absolute deviances).
    // void filterByPeakWidth(std::vector<MassTrace>&, std::vector<MassTrace>&);

    /** @name Main computation methods
        */
    /// Main method of MassTraceDetection. Extracts mass traces of a @ref MSExperiment and gathers them into a vector container.
    void run(const MSExperiment<Peak1D> &, std::vector<MassTrace> &);

    /// Invokes the run method (see above) on merely a subregion of a @ref MSExperiment map.
    void run(MSExperiment<Peak1D>::ConstAreaIterator & begin, MSExperiment<Peak1D>::ConstAreaIterator & end, std::vector<MassTrace> & found_masstraces);

protected:
    virtual void updateMembers_();

private:
    // parameter stuff
    DoubleReal mass_error_ppm_;
    DoubleReal noise_threshold_int_;
    DoubleReal chrom_peak_snr_;

    DoubleReal min_sample_rate_;
    DoubleReal min_peak_width_;

    bool reestimate_mt_sd_;
  };
}

#endif // OPENMS_FILTERING_DATAREDUCTION_MASSTRACEDETECTION_H
