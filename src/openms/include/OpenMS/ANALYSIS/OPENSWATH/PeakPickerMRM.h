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

#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

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

    //@{
    /// Constructor
    PeakPickerMRM();

    /// Destructor
    ~PeakPickerMRM() override {}
    //@}

    /**
      @brief Finds peaks in a single chromatogram and annotates left/right borders

      It uses a modified algorithm of the PeakPickerHiRes

      This function will return a picked chromatogram
    */
    void pickChromatogram(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom);
    

    /**
      @brief Finds peaks in a single chromatogram and annotates left/right borders

      It uses a modified algorithm of the PeakPickerHiRes

      This function will return a picked chromatogram and a smoothed chromatogram
    */
    void pickChromatogram(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom, MSChromatogram& smoothed_chrom);

protected:

    void pickChromatogramCrawdad_(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom);

    void pickChromatogram_(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom);

    /**
      @brief Compute peak area (peak integration)
    */
    void integratePeaks_(const MSChromatogram& chromatogram);

    /**
      @brief Helper function to find the closest peak in a chromatogram to "target_rt" 

      The search will start from the index current_peak, so the function is
      assuming the closest peak is to the right of current_peak.

      It will return the index of the closest peak in the chromatogram.
    */
    Size findClosestPeak_(const MSChromatogram& chromatogram, double target_rt, Size current_peak = 0);

    /**
      @brief Helper function to remove overlapping peaks in a single Chromatogram

    */
    void removeOverlappingPeaks_(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom);


    /// Synchronize members with param class
    void updateMembers_() override;

    /// Assignment operator is protected for algorithm
    PeakPickerMRM& operator=(const PeakPickerMRM& rhs);

    // Members
    /// Frame length for the SGolay smoothing
    UInt sgolay_frame_length_;
    /// Polynomial order for the SGolay smoothing
    UInt sgolay_polynomial_order_;
    /// Width of the Gaussian smoothing
    double gauss_width_;
    /// Whether to use Gaussian smoothing
    bool use_gauss_;
    /// Whether to resolve overlapping peaks 
    bool remove_overlapping_;

    /// Forced peak with
    double peak_width_;
    /// Signal to noise threshold
    double signal_to_noise_;

    /// Signal to noise window length
    double sn_win_len_;
    /// Signal to noise bin count
    UInt sn_bin_count_;
    /// Whether to write out log messages of the SN estimator
    bool write_sn_log_messages_;
    /// Peak picker method
    String method_;

    /// Temporary vector to hold the integrated intensities
    std::vector<double> integrated_intensities_;
    /// Temporary vector to hold the peak left widths
    std::vector<int> left_width_;
    /// Temporary vector to hold the peak right widths
    std::vector<int> right_width_;

    PeakPickerHiRes pp_;
    SavitzkyGolayFilter sgolay_;
    GaussFilter gauss_;
    SignalToNoiseEstimatorMedian<MSChromatogram > snt_;
  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_PEAKPICKERMRM_H

