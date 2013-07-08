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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_PEAKPICKERMRM_H
#define OPENMS_ANALYSIS_OPENSWATH_PEAKPICKERMRM_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>

#ifdef WITH_CRAWDAD
#include <CrawdadWrapper.h>
#endif

namespace OpenMS
{

  /**
  @brief The PeakPickerMRM finds peaks a single chromatogram.

  @htmlinclude OpenMS_PeakPickerMRM.parameters

  It uses the PeakPickerHiRes internally to find interesting seed candidates.
  These candidates are then expanded and a right/left border of the peak is
  searched.
  Additionally, overlapping peaks can be removed.

  */

  class OPENMS_DLLAPI PeakPickerMRM :
    public DefaultParamHandler
  {

public:

    // this is the type in which we store the chromatograms for this analysis
    typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram; 

    //@{
    /// Constructor
    PeakPickerMRM();

    /// Destructor
    ~PeakPickerMRM() {}
    //@}

    /**
      @brief Finds peaks in a single chromatogram and annotates left/right borders

      It uses a modified algorithm of the PeakPickerHiRes

      This function will return a smoothed chromatogram and a picked chromatogram
    */
    void pickChromatogram(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& picked_chrom);

protected:

    void pickChromatogramCrowdad(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& picked_chrom);

    void pickChromatogram_(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& picked_chrom);

    /**
      @brief Compute peak area (peak integration)
    */
    void integratePeaks_(const RichPeakChromatogram& chromatogram);

    /**
      @brief Helper function to find the closest peak in a chromatogram to "target_rt" 

      The search will start from the index current_peak, so the function is
      assuming the closest peak is to the right of current_peak.

      It will return the index of the closest peak in the chromatogram.
    */
    Size findClosestPeak_(const RichPeakChromatogram& chromatogram, double target_rt, Size current_peak = 0);

    /**
      @brief Helper function to remove overlapping peaks in a single Chromatogram

    */
    void removeOverlappingPeaks_(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& picked_chrom);


    /// Synchronize members with param class
    void updateMembers_();

    /// Assignment operator is protected for algorithm
    PeakPickerMRM& operator=(const PeakPickerMRM& rhs);

    // Members
    UInt sgolay_frame_length_;
    UInt sgolay_polynomial_order_;
    DoubleReal gauss_width_;
    bool use_gauss_;
    bool remove_overlapping_;

    DoubleReal peak_width_;
    DoubleReal signal_to_noise_;

    DoubleReal sn_win_len_;
    UInt sn_bin_count_;
    String method_;

    std::vector<double> integrated_intensities;
    std::vector<int> left_width;
    std::vector<int> right_width;

  };
}

#endif
